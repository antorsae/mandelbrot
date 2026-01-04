// test_flags.cpp - Test rendering correctness across optimization flag combinations
// Compares optimized rendering against fully-disabled baseline
// Uses --dump-iter to get raw iteration buffers for comparison
//
// Enhanced per Codex recommendations:
// - Unique temp files per run (avoids stale data on crashes)
// - Exit code checking
// - Higher resolution option (--highres)
// - Cache-hit test

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <unistd.h>
#include <sys/wait.h>

// Global test resolution (can be overridden with --highres)
static int g_test_width = 64;
static int g_test_height = 40;
static bool g_verbose = false;

// Test location
struct TestLocation {
    const char* name;
    const char* pos;       // Complex position string
    double zoom;
    int expected_max_iter; // Approximate, for validation
};

// Test results
struct TestResult {
    std::string location;
    std::string flags;
    int width, height, max_iter;
    int total_pixels;
    int escape_match;      // Escaped/bounded classification matches
    int exact_match;       // Iteration values match exactly
    double max_error;      // Max |baseline - optimized| for escaped points
    double mean_error;     // Mean |baseline - optimized| for escaped points
    int block_artifacts;   // 16x16 blocks with >50% mismatch (tile artifacts)
    bool passed;
};

// Generate unique temp file path
std::string make_temp_path() {
    char path[256];
    snprintf(path, sizeof(path), "/tmp/test_flags_%d_%ld.bin", getpid(), random());
    return path;
}

// Read iteration buffer from file
std::vector<double> read_iter_buffer(const char* path, int& width, int& height, int& max_iter) {
    FILE* fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "ERROR: Cannot open %s\n", path);
        return {};
    }

    int32_t header[3];
    if (fread(header, sizeof(int32_t), 3, fp) != 3) {
        fprintf(stderr, "ERROR: Cannot read header from %s\n", path);
        fclose(fp);
        return {};
    }

    width = header[0];
    height = header[1];
    max_iter = header[2];

    std::vector<double> data(width * height);
    size_t read = fread(data.data(), sizeof(double), data.size(), fp);
    fclose(fp);

    if (read != data.size()) {
        fprintf(stderr, "ERROR: Short read from %s (got %zu, expected %zu)\n",
                path, read, data.size());
        return {};
    }

    return data;
}

// Run mandelbrot with given flags and return iteration buffer
// Returns empty vector on failure
std::vector<double> run_mandelbrot(const TestLocation& loc, const char* flags,
                                    int& width, int& height, int& max_iter,
                                    int req_width = 0, int req_height = 0) {
    if (req_width == 0) req_width = g_test_width;
    if (req_height == 0) req_height = g_test_height;

    std::string tmpfile = make_temp_path();
    char cmd[4096];

    // Use explicit dimensions for consistent sizing
    snprintf(cmd, sizeof(cmd),
             "./mandelbrot --pos '%s' --zoom %.17g %s --width %d --height %d --dump-iter=%s 2>&1",
             loc.pos, loc.zoom, flags, req_width, req_height, tmpfile.c_str());

    if (g_verbose) {
        printf("    CMD: %s\n", cmd);
    }

    FILE* pipe = popen(cmd, "r");
    if (!pipe) {
        fprintf(stderr, "ERROR: Cannot run mandelbrot\n");
        return {};
    }

    char output[1024];
    while (fgets(output, sizeof(output), pipe)) {
        if (g_verbose) {
            printf("    OUT: %s", output);
        }
    }

    int status = pclose(pipe);
    int exit_code = WEXITSTATUS(status);

    if (exit_code != 0) {
        fprintf(stderr, "ERROR: mandelbrot exited with code %d\n", exit_code);
        unlink(tmpfile.c_str());
        return {};
    }

    auto result = read_iter_buffer(tmpfile.c_str(), width, height, max_iter);

    // Clean up temp file
    unlink(tmpfile.c_str());

    return result;
}

// Compare two iteration buffers
TestResult compare_buffers(const std::vector<double>& baseline,
                           const std::vector<double>& optimized,
                           int width, int height, int max_iter,
                           const char* location, const char* flags) {
    TestResult result;
    result.location = location;
    result.flags = flags;
    result.width = width;
    result.height = height;
    result.max_iter = max_iter;
    result.total_pixels = width * height;
    result.escape_match = 0;
    result.exact_match = 0;
    result.max_error = 0.0;
    result.mean_error = 0.0;
    result.block_artifacts = 0;
    result.passed = false;

    if (baseline.size() != optimized.size()) {
        fprintf(stderr, "ERROR: Buffer size mismatch\n");
        return result;
    }

    int escaped_count = 0;
    double error_sum = 0.0;

    for (size_t i = 0; i < baseline.size(); i++) {
        double b = baseline[i];
        double o = optimized[i];

        // -1.0 means bounded (interior), >= 0 means escaped with smooth iteration
        bool b_bounded = (b < 0);
        bool o_bounded = (o < 0);

        // Escape/bounded classification match?
        if (b_bounded == o_bounded) {
            result.escape_match++;

            // For escaped points, check iteration value
            if (!b_bounded) {
                escaped_count++;
                double err = std::abs(b - o);
                if (err < 0.001) {
                    result.exact_match++;
                }
                error_sum += err;
                if (err > result.max_error) {
                    result.max_error = err;
                }
            } else {
                // Both bounded - exact match
                result.exact_match++;
            }
        }
    }

    if (escaped_count > 0) {
        result.mean_error = error_sum / escaped_count;
    }

    // Check for block artifacts (16x16 tiles with many mismatches)
    const int TILE_SIZE = 16;
    for (int ty = 0; ty < height; ty += TILE_SIZE) {
        for (int tx = 0; tx < width; tx += TILE_SIZE) {
            int tile_total = 0;
            int tile_mismatch = 0;

            for (int y = ty; y < ty + TILE_SIZE && y < height; y++) {
                for (int x = tx; x < tx + TILE_SIZE && x < width; x++) {
                    tile_total++;
                    size_t idx = y * width + x;
                    double b = baseline[idx];
                    double o = optimized[idx];

                    bool b_bounded = (b < 0);
                    bool o_bounded = (o < 0);

                    if (b_bounded != o_bounded ||
                        (!b_bounded && std::abs(b - o) > 1.0)) {
                        tile_mismatch++;
                    }
                }
            }

            // Flag tile as artifact if >50% pixels have significant differences
            if (tile_mismatch > tile_total / 2) {
                result.block_artifacts++;
            }
        }
    }

    // Pass criteria:
    // - At least 99% escape/bounded classification match
    // - Max error < 2.0 for smooth iteration values
    // - No block artifacts
    double escape_pct = 100.0 * result.escape_match / result.total_pixels;
    result.passed = (escape_pct >= 99.0) &&
                    (result.max_error < 2.0) &&
                    (result.block_artifacts == 0);

    return result;
}

// Print test result
void print_result(const TestResult& r) {
    double escape_pct = 100.0 * r.escape_match / r.total_pixels;
    double exact_pct = 100.0 * r.exact_match / r.total_pixels;

    printf("  %-20s %-35s %s\n", r.location.c_str(), r.flags.c_str(),
           r.passed ? "PASS" : "FAIL");
    printf("    Escape match: %d/%d (%.1f%%)  Exact: %.1f%%  MaxErr: %.2f  MeanErr: %.3f  BlockArt: %d\n",
           r.escape_match, r.total_pixels, escape_pct, exact_pct,
           r.max_error, r.mean_error, r.block_artifacts);
}

// Cache-hit test: render same frame twice with cache enabled, compare results
// This tests that the orbit cache produces identical results
bool test_cache_hit(const TestLocation& loc) {
    printf("\n[Cache-hit test: %s]\n", loc.name);

    // First render (cache miss - fresh orbit computation)
    int w1, h1, max1;
    auto buf1 = run_mandelbrot(loc, "", w1, h1, max1);
    if (buf1.empty()) {
        printf("  FAIL: First render failed\n");
        return false;
    }

    // Second render (should be cache hit if same position)
    // Note: Each process is independent, so cache hit won't happen across processes.
    // This test validates that the --no-cache flag produces same results as cached.
    int w2, h2, max2;
    auto buf2 = run_mandelbrot(loc, "--no-cache", w2, h2, max2);
    if (buf2.empty()) {
        printf("  FAIL: Second render (--no-cache) failed\n");
        return false;
    }

    // Compare
    auto result = compare_buffers(buf1, buf2, w1, h1, max1, loc.name, "cache vs no-cache");

    double exact_pct = 100.0 * result.exact_match / result.total_pixels;
    if (exact_pct == 100.0) {
        printf("  PASS: Cache vs no-cache identical (100%% exact match)\n");
        return true;
    } else {
        printf("  FAIL: Cache vs no-cache differ (%.1f%% exact, max_err=%.3f)\n",
               exact_pct, result.max_error);
        return false;
    }
}

void print_usage(const char* prog) {
    printf("Usage: %s [options]\n", prog);
    printf("Options:\n");
    printf("  --highres    Use 256x160 resolution instead of 64x40\n");
    printf("  --fullres    Use 640x400 resolution (slow but matches --image default)\n");
    printf("  --verbose    Print commands and output\n");
    printf("  --cache-only Only run cache-hit tests\n");
    printf("  --help       Show this help\n");
}

int main(int argc, char** argv) {
    bool cache_only = false;

    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--highres") == 0) {
            g_test_width = 256;
            g_test_height = 160;
        } else if (strcmp(argv[i], "--fullres") == 0) {
            g_test_width = 640;
            g_test_height = 400;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            g_verbose = true;
        } else if (strcmp(argv[i], "--cache-only") == 0) {
            cache_only = true;
        } else if (strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    // Seed random for unique temp files
    srandom(time(nullptr) ^ getpid());

    printf("=== Mandelbrot Optimization Flag Test Suite ===\n");
    printf("Resolution: %dx%d\n\n", g_test_width, g_test_height);

    // Test locations per Codex recommendations
    std::vector<TestLocation> locations = {
        // Seahorse Valley - classic test location
        {"Seahorse-1e6", "-0.743643887037151+0.131825904205330i", 1e6, 500},
        {"Seahorse-1e10", "-0.743643887037151+0.131825904205330i", 1e10, 700},
        {"Seahorse-1e12", "-0.743643887037151+0.131825904205330i", 1e12, 850},

        // Original bug report location
        {"BugReport-3e12", "-0.7307997873-0.1763509361i", 3.18e12, 900},

        // Boundary minibrot region
        {"Minibrot-2e12", "-0.7210622206-0.1999961551i", 2.12e12, 900},

        // Uniform interior (should be all bounded)
        {"Interior-1e6", "0.0+0.0i", 1e6, 500},

        // Far exterior (should be all escaped quickly)
        {"Exterior-1e6", "0.5+0.0i", 1e6, 500},

        // The extreme zoom artifact location (from user bug report)
        {"Extreme-1e18", "-0.99999990959787322:7.6610371875867799e-18-0.29019415957242539:-1.5207291758725448e-17i", 1.2324511097589686e18, 1500},
    };

    // Cache-hit tests
    printf("=== Cache Consistency Tests ===\n");
    int cache_pass = 0, cache_fail = 0;
    for (const auto& loc : locations) {
        if (test_cache_hit(loc)) {
            cache_pass++;
        } else {
            cache_fail++;
        }
    }
    printf("\nCache tests: %d passed, %d failed\n", cache_pass, cache_fail);

    if (cache_only) {
        return cache_fail > 0 ? 1 : 0;
    }

    printf("\n=== Flag Combination Tests ===\n");

    // Flag combinations to test
    // Baseline: all optimizations disabled
    // Then test each optimization enabled individually and all enabled
    std::vector<std::pair<const char*, const char*>> flag_sets = {
        {"all-off", "--no-sa --no-bla --no-block --no-cache"},  // Baseline
        {"sa-only", "--no-bla --no-block --no-cache"},          // Only SA enabled
        {"bla-only", "--no-sa --no-block --no-cache"},          // Only BLA enabled (now works!)
        {"block-only", "--no-sa --no-bla --no-cache"},          // Only block enabled
        {"cache-only", "--no-sa --no-bla --no-block"},          // Only cache enabled
        {"sa+bla", "--no-block --no-cache"},                    // SA + BLA
        {"bla+block", "--no-sa --no-cache"},                    // BLA + block (no SA)
        {"sa+block", "--no-bla --no-cache"},                    // SA + block
        {"all-on", ""},                                          // All optimizations
    };

    std::vector<TestResult> results;
    int total_pass = 0;
    int total_fail = 0;

    for (const auto& loc : locations) {
        printf("\nTesting: %s (zoom=%.2e)\n", loc.name, loc.zoom);
        printf("─────────────────────────────────────────────────────────────────────\n");

        // Get baseline (all optimizations off)
        int base_w, base_h, base_max;
        auto baseline = run_mandelbrot(loc, flag_sets[0].second, base_w, base_h, base_max);

        if (baseline.empty()) {
            printf("  ERROR: Failed to get baseline for %s\n", loc.name);
            total_fail++;
            continue;
        }

        printf("  Baseline: %dx%d, max_iter=%d\n", base_w, base_h, base_max);

        // Test each flag combination against baseline
        for (size_t i = 1; i < flag_sets.size(); i++) {
            int opt_w, opt_h, opt_max;
            auto optimized = run_mandelbrot(loc, flag_sets[i].second, opt_w, opt_h, opt_max);

            if (optimized.empty()) {
                printf("  ERROR: Failed to get result for %s with %s\n",
                       loc.name, flag_sets[i].first);
                total_fail++;
                continue;
            }

            auto result = compare_buffers(baseline, optimized, base_w, base_h, base_max,
                                          loc.name, flag_sets[i].first);
            print_result(result);
            results.push_back(result);

            if (result.passed) {
                total_pass++;
            } else {
                total_fail++;
            }
        }
    }

    // Summary
    printf("\n═══════════════════════════════════════════════════════════════════════\n");
    printf("SUMMARY:\n");
    printf("  Cache tests: %d passed, %d failed\n", cache_pass, cache_fail);
    printf("  Flag tests:  %d passed, %d failed\n", total_pass, total_fail);
    printf("  Total:       %d passed, %d failed\n",
           cache_pass + total_pass, cache_fail + total_fail);
    printf("═══════════════════════════════════════════════════════════════════════\n");

    // List failures
    int all_fail = cache_fail + total_fail;
    if (all_fail > 0) {
        printf("\nFailed tests:\n");
        for (const auto& r : results) {
            if (!r.passed) {
                printf("  - %s @ %s\n", r.location.c_str(), r.flags.c_str());
            }
        }
    }

    return all_fail > 0 ? 1 : 0;
}
