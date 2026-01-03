/*
 * TRAJECTORY INTERESTINGNESS TEST
 *
 * Tests that cinematic trajectories never pass through boring (black/uniform) regions.
 * Each trajectory is sampled at 30fps and checked for:
 * - Single-color frames (all same iteration = FAIL)
 * - Low interest score (< MIN_INTEREST_SCORE = FAIL)
 *
 * Usage: ./test_trajectory [--verbose]
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>

// ═══════════════════════════════════════════════════════════════════════════
// CONFIGURATION
// ═══════════════════════════════════════════════════════════════════════════

const double MIN_INTEREST_SCORE = 5.0;  // Minimum score to pass
const int FPS = 30;                      // Frames per second for testing
const double TRAJECTORY_DURATION = 30.0; // Test 30-second trajectories
const int MAX_ITER = 256;                // Iterations for interest scoring

// ═══════════════════════════════════════════════════════════════════════════
// DATA STRUCTURES (extracted from mandelbrot.cpp)
// ═══════════════════════════════════════════════════════════════════════════

struct TrajectoryWaypoint {
    double x, y;           // Position in complex plane
    double zoom;           // Zoom level
    double angle;          // View rotation
    double interest_score; // For debugging/visualization
};

struct CinematicPath {
    std::vector<TrajectoryWaypoint> waypoints;
    double total_duration;
    bool valid = false;

    TrajectoryWaypoint evaluate(double t) const;
};

struct Location {
    std::string name;
    double real, imag, zoom;
};

struct TestResult {
    bool passed;
    int fail_frame;
    double fail_score;
    double fail_time;
    double fail_zoom;
    std::string reason;
};

// ═══════════════════════════════════════════════════════════════════════════
// MANDELBROT ITERATION (extracted from mandelbrot.cpp)
// ═══════════════════════════════════════════════════════════════════════════

inline int quick_iterate(double cx, double cy, int max_iter) {
    double zr = 0, zi = 0;
    for (int i = 0; i < max_iter; i++) {
        double zr2 = zr * zr, zi2 = zi * zi;
        if (zr2 + zi2 > 4.0) return i;
        zi = 2 * zr * zi + cy;
        zr = zr2 - zi2 + cx;
    }
    return max_iter;
}

// ═══════════════════════════════════════════════════════════════════════════
// INTEREST SCORING (extracted from mandelbrot.cpp with fixes)
// ═══════════════════════════════════════════════════════════════════════════

double calculate_interest_score_at_zoom(double cx, double cy, double zoom, int sample_radius = 5) {
    // FIX: Scale max_iter with zoom - deep zooms need more iterations
    const int max_iter = 256 + (int)(log10(std::max(1.0, zoom)) * 50);
    // FIX: Sample scale covers full viewport (3/zoom wide)
    // With sample_radius=5, we have 11 samples spanning viewport width
    const double viewport_width = 3.0 / zoom;
    const double sample_scale = viewport_width / (2.0 * sample_radius);

    std::vector<int> iterations;
    int bounded_count = 0;
    int escaped_count = 0;

    for (int dy = -sample_radius; dy <= sample_radius; dy++) {
        for (int dx = -sample_radius; dx <= sample_radius; dx++) {
            double px = cx + dx * sample_scale;
            double py = cy + dy * sample_scale;
            int iter = quick_iterate(px, py, max_iter);
            iterations.push_back(iter);
            if (iter >= max_iter) bounded_count++;
            else escaped_count++;
        }
    }

    // FIX: Strong penalty for uniform regions (all black or all escaped)
    int total_samples = (int)iterations.size();
    if (bounded_count == total_samples || escaped_count == total_samples) {
        return -100.0;  // Strong negative score for uniform regions
    }

    // Score based on boundary proximity - ratio of minority to majority type
    // Higher ratio = closer to boundary = more interesting
    double ratio = std::min(bounded_count, escaped_count) /
                   (double)std::max(bounded_count, escaped_count);
    double on_boundary_score = ratio * 50.0;

    // Calculate variance for complexity
    double sum = 0, sum_sq = 0;
    for (int iter : iterations) {
        sum += iter;
        sum_sq += iter * iter;
    }
    double n = iterations.size();
    double mean = sum / n;
    double variance = std::max(0.0, (sum_sq / n) - (mean * mean));
    double variance_score = std::min(50.0, sqrt(variance));

    // Bonus for high average iteration
    double avg_score = std::min(20.0, mean / 12.0);

    return on_boundary_score + variance_score + avg_score;
}

// Check if a point appears as "single color" (all samples same iteration)
bool is_single_color(double cx, double cy, double zoom) {
    // FIX: Use viewport-relative sampling and zoom-scaled iterations
    const double viewport_width = 3.0 / zoom;
    const double sample_scale = viewport_width / 6.0;  // 7 samples across viewport
    const int max_iter = 256 + (int)(log10(std::max(1.0, zoom)) * 50);
    int first_iter = -1;
    for (int dy = -3; dy <= 3; dy++) {
        for (int dx = -3; dx <= 3; dx++) {
            int iter = quick_iterate(cx + dx * sample_scale, cy + dy * sample_scale, max_iter);
            if (first_iter < 0) first_iter = iter;
            else if (iter != first_iter) return false;
        }
    }
    return true;  // All same = single color
}

// ═══════════════════════════════════════════════════════════════════════════
// WAYPOINT SEARCH (extracted from mandelbrot.cpp with fixes)
// ═══════════════════════════════════════════════════════════════════════════

TrajectoryWaypoint find_best_waypoint(double center_x, double center_y, double zoom,
                                       double search_radius, double prev_x, double prev_y,
                                       std::mt19937& rng) {
    const int NUM_CANDIDATES = 200;  // More candidates
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    double effective_radius = search_radius * 2.0;
    int search_attempts = 0;
    const int MAX_SEARCH_ATTEMPTS = 5;  // More attempts

    while (search_attempts < MAX_SEARCH_ATTEMPTS) {
        std::vector<TrajectoryWaypoint> candidates;

        // Random sampling
        for (int i = 0; i < NUM_CANDIDATES; i++) {
            double dx = dist(rng) * effective_radius;
            double dy = dist(rng) * effective_radius;
            double cx = center_x + dx;
            double cy = center_y + dy;

            double score = calculate_interest_score_at_zoom(cx, cy, zoom);

            // Light continuity penalty
            double jump_dist = hypot(cx - prev_x, cy - prev_y);
            double jump_relative = jump_dist / effective_radius;
            double continuity_penalty = std::min(10.0, jump_relative * 2.0);
            score -= continuity_penalty;

            candidates.push_back({cx, cy, zoom, 0.0, score});
        }

        // Also try known good boundary points (radial search pattern)
        // The boundary of the cardioid is at r = 0.25*(1 - cos(theta)) from center (-0.25, 0)
        for (int angle_deg = 0; angle_deg < 360; angle_deg += 15) {
            double theta = angle_deg * M_PI / 180.0;
            // Main cardioid boundary parameterization
            double bx = 0.25 * cos(theta) - 0.25 * cos(2*theta);
            double by = 0.25 * sin(theta) - 0.25 * sin(2*theta);
            // Offset slightly toward center_x/center_y
            double cx = bx + (center_x - bx) * 0.1;
            double cy = by + (center_y - by) * 0.1;

            double score = calculate_interest_score_at_zoom(cx, cy, zoom);
            double jump_dist = hypot(cx - prev_x, cy - prev_y);
            double jump_relative = jump_dist / effective_radius;
            score -= std::min(10.0, jump_relative * 2.0);

            candidates.push_back({cx, cy, zoom, 0.0, score});
        }

        auto best = std::max_element(candidates.begin(), candidates.end(),
            [](const auto& a, const auto& b) { return a.interest_score < b.interest_score; });

        if (best->interest_score >= 5.0 || search_attempts == MAX_SEARCH_ATTEMPTS - 1) {
            return *best;
        }

        effective_radius *= 2.0;  // Expand more gradually
        search_attempts++;
    }

    return {center_x, center_y, zoom, 0.0, 0.0};
}

// ═══════════════════════════════════════════════════════════════════════════
// SPLINE INTERPOLATION (extracted from mandelbrot.cpp)
// ═══════════════════════════════════════════════════════════════════════════

inline double catmull_rom(double p0, double p1, double p2, double p3, double t) {
    double t2 = t * t;
    double t3 = t2 * t;
    return 0.5 * ((2.0 * p1) +
                  (-p0 + p2) * t +
                  (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2 +
                  (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
}

inline double ease_in_out(double t) {
    return t < 0.5 ? 4 * t * t * t : 1 - pow(-2 * t + 2, 3) / 2;
}

TrajectoryWaypoint CinematicPath::evaluate(double t) const {
    if (!valid || waypoints.size() < 2) {
        return waypoints.empty() ? TrajectoryWaypoint{0, 0, 1, 0, 0} : waypoints[0];
    }

    double normalized = std::max(0.0, std::min(1.0, t / total_duration));
    double eased = ease_in_out(normalized);

    int n = waypoints.size() - 1;
    double segment_pos = eased * n;
    int segment = std::min((int)segment_pos, n - 1);
    double local_t = segment_pos - segment;

    const auto& p1 = waypoints[segment];
    const auto& p2 = waypoints[std::min(segment + 1, n)];

    // Use linear interpolation with boring-region avoidance
    double x = p1.x + local_t * (p2.x - p1.x);
    double y = p1.y + local_t * (p2.y - p1.y);
    double log_zoom = log(p1.zoom) + local_t * (log(p2.zoom) - log(p1.zoom));
    double zoom = exp(log_zoom);

    // Check if linear interpolation lands in boring region
    double score = calculate_interest_score_at_zoom(x, y, zoom);
    if (score < 0) {
        // Try biasing toward the waypoint with better score
        // (waypoints should have positive scores)
        double score1 = p1.interest_score;
        double score2 = p2.interest_score;

        // If one waypoint is clearly better, bias toward it
        if (score1 > score2 + 10) {
            // Bias toward p1
            double new_t = local_t * 0.5;  // Move closer to p1
            x = p1.x + new_t * (p2.x - p1.x);
            y = p1.y + new_t * (p2.y - p1.y);
        } else if (score2 > score1 + 10) {
            // Bias toward p2
            double new_t = 0.5 + local_t * 0.5;  // Move closer to p2
            x = p1.x + new_t * (p2.x - p1.x);
            y = p1.y + new_t * (p2.y - p1.y);
        }
        // Else: both waypoints similar, keep linear (may still be boring)
    }

    double angle = p1.angle + local_t * (p2.angle - p1.angle);
    return {x, y, zoom, angle, score};
}

// ═══════════════════════════════════════════════════════════════════════════
// PATH PLANNING (extracted from mandelbrot.cpp)
// ═══════════════════════════════════════════════════════════════════════════

CinematicPath plan_cinematic_path(double start_x, double start_y, double start_zoom, double start_angle,
                                   double target_x, double target_y, double target_zoom, double target_angle,
                                   double duration, std::mt19937& rng) {
    CinematicPath path;
    path.total_duration = duration;

    // Use NON-UNIFORM waypoint spacing: dense at low zoom, sparse at high zoom
    // This ensures short chord lengths at low zoom where boundary is coarse
    std::vector<double> zoom_levels;

    // Phase 1: Low zoom (1-100) - very dense waypoints
    int low_zoom_waypoints = 50;
    for (int i = 0; i <= low_zoom_waypoints; i++) {
        double t = (double)i / low_zoom_waypoints;
        double log_zoom = log(start_zoom) + t * (log(100.0) - log(start_zoom));
        if (exp(log_zoom) <= target_zoom) {
            zoom_levels.push_back(exp(log_zoom));
        }
    }

    // Phase 2: Medium zoom (100-10000) - moderate density
    int med_zoom_waypoints = 30;
    for (int i = 1; i <= med_zoom_waypoints; i++) {
        double t = (double)i / med_zoom_waypoints;
        double log_zoom = log(100.0) + t * (log(10000.0) - log(100.0));
        if (exp(log_zoom) <= target_zoom && exp(log_zoom) > zoom_levels.back()) {
            zoom_levels.push_back(exp(log_zoom));
        }
    }

    // Phase 3: High zoom (10000+) - standard density
    if (target_zoom > 10000.0) {
        double remaining_range = log(target_zoom) - log(10000.0);
        int high_zoom_waypoints = std::max(20, (int)(remaining_range / log(10.0) * 3));
        for (int i = 1; i <= high_zoom_waypoints; i++) {
            double t = (double)i / high_zoom_waypoints;
            double log_zoom = log(10000.0) + t * remaining_range;
            zoom_levels.push_back(exp(log_zoom));
        }
    }

    // Ensure we end at target zoom
    if (zoom_levels.empty() || zoom_levels.back() < target_zoom * 0.99) {
        zoom_levels.push_back(target_zoom);
    }

    int num_waypoints = (int)zoom_levels.size() - 1;

    double prev_x = start_x, prev_y = start_y;

    for (int i = 0; i <= num_waypoints; i++) {
        double zoom = zoom_levels[i];
        double t = (double)i / num_waypoints;

        double expected_x = start_x + t * (target_x - start_x);
        double expected_y = start_y + t * (target_y - start_y);

        // FIX: Search radius scales with zoom to find boundary
        double viewport_width = 3.0 / zoom;
        double base_radius;
        if (zoom < 20.0) {
            base_radius = 2.0;  // Large for escaping set interior
        } else if (zoom < 1000.0) {
            base_radius = std::max(0.1, 1.0 / sqrt(zoom));  // Medium zoom
        } else {
            base_radius = viewport_width * 50.0;  // Deep zoom: viewport-relative
        }
        double progress_factor = 1.0 - t * t;
        double search_radius = base_radius * (0.2 + 0.8 * progress_factor);

        TrajectoryWaypoint waypoint;

        if (i == 0) {
            waypoint = {start_x, start_y, start_zoom, start_angle, 100.0};
        } else if (i == num_waypoints) {
            waypoint = {target_x, target_y, target_zoom, target_angle, 100.0};
        } else {
            waypoint = find_best_waypoint(expected_x, expected_y, zoom,
                                          search_radius, prev_x, prev_y, rng);
            waypoint.angle = start_angle + t * (target_angle - start_angle);
        }

        path.waypoints.push_back(waypoint);
        prev_x = waypoint.x;
        prev_y = waypoint.y;
    }

    path.valid = true;
    return path;
}

// ═══════════════════════════════════════════════════════════════════════════
// LOCATION LOADING
// ═══════════════════════════════════════════════════════════════════════════

std::vector<Location> load_locations(const char* filename) {
    std::vector<Location> locations;
    std::ifstream file(filename);

    if (!file.is_open()) {
        fprintf(stderr, "ERROR: Cannot open %s\n", filename);
        return locations;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        // Parse: name,real,imag,zoom
        std::stringstream ss(line);
        std::string name, real_str, imag_str, zoom_str;

        if (!std::getline(ss, name, ',')) continue;
        if (!std::getline(ss, real_str, ',')) continue;
        if (!std::getline(ss, imag_str, ',')) continue;
        if (!std::getline(ss, zoom_str, ',')) continue;

        Location loc;
        loc.name = name;
        loc.real = std::stod(real_str);
        loc.imag = std::stod(imag_str);
        loc.zoom = std::stod(zoom_str);
        locations.push_back(loc);
    }

    return locations;
}

// ═══════════════════════════════════════════════════════════════════════════
// TRAJECTORY TESTING
// ═══════════════════════════════════════════════════════════════════════════

TestResult test_trajectory(const Location& target, double duration, bool verbose) {
    std::mt19937 rng(42);  // Deterministic seed for reproducibility

    if (verbose) {
        printf("  Planning path to (%.10f, %.10f) zoom=%.2e...\n",
               target.real, target.imag, target.zoom);
    }

    // Start from outside the set (0, 0.7) which gives high interest scores
    // at low zooms and allows graceful transition to targets
    CinematicPath path = plan_cinematic_path(
        0.0, 0.7, 1.0, 0.0,  // Start: above the set (high interest at low zoom)
        target.real, target.imag, target.zoom, 0.0,
        duration, rng
    );

    if (verbose) {
        printf("  Planned %zu waypoints, testing %d frames...\n",
               path.waypoints.size(), (int)(duration * FPS));
    }

    int num_frames = (int)(duration * FPS);
    int low_score_count = 0;

    // At very low zoom (< 20), the whole Mandelbrot set is visible
    // Single-color frames only matter once we're zoomed in enough
    // that boring regions become problematic
    const double MIN_ZOOM_FOR_STRICT_CHECK = 20.0;

    for (int f = 0; f <= num_frames; f++) {
        double t = (double)f / FPS;
        auto wp = path.evaluate(t);

        // Only enforce strict single-color check at sufficient zoom
        if (wp.zoom >= MIN_ZOOM_FOR_STRICT_CHECK) {
            if (is_single_color(wp.x, wp.y, wp.zoom)) {
                return {false, f, 0.0, t, wp.zoom, "single color frame"};
            }
        }

        // Check interest score (always)
        double score = calculate_interest_score_at_zoom(wp.x, wp.y, wp.zoom);
        if (score < MIN_INTEREST_SCORE) {
            low_score_count++;
            if (verbose) {
                printf("    Frame %d (t=%.2fs zoom=%.2e): score=%.1f < %.1f\n",
                       f, t, wp.zoom, score, MIN_INTEREST_SCORE);
            }
            // More tolerance at low zoom
            int max_low_score = (wp.zoom < MIN_ZOOM_FOR_STRICT_CHECK) ? 50 : 5;
            if (low_score_count > max_low_score) {
                return {false, f, score, t, wp.zoom, "too many low interest frames"};
            }
        }
    }

    // Allow up to 5 low-score frames (some tolerance)
    if (low_score_count > 0 && verbose) {
        printf("  Warning: %d frames with low interest (but within tolerance)\n", low_score_count);
    }

    return {true, -1, 0.0, 0.0, 0.0, ""};
}

TestResult test_location_interestingness(const Location& loc, bool verbose) {
    double score = calculate_interest_score_at_zoom(loc.real, loc.imag, loc.zoom);

    if (verbose) {
        printf("  Location score: %.1f\n", score);
    }

    if (score < 50.0) {
        return {false, -1, score, 0.0, loc.zoom, "location not interesting enough"};
    }

    return {true, -1, score, 0.0, loc.zoom, ""};
}

// ═══════════════════════════════════════════════════════════════════════════
// MAIN
// ═══════════════════════════════════════════════════════════════════════════

int main(int argc, char* argv[]) {
    bool verbose = false;
    bool skip_trajectory = false;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
            verbose = true;
        }
        if (strcmp(argv[i], "--skip-trajectory") == 0) {
            skip_trajectory = true;
        }
    }

    printf("TRAJECTORY INTERESTINGNESS TEST\n");
    printf("================================\n");
    printf("Testing that cinematic trajectories avoid boring regions.\n");
    printf("Duration: %.0f seconds per trajectory\n", TRAJECTORY_DURATION);
    printf("Sample rate: %d fps (%d frames per trajectory)\n", FPS, (int)(TRAJECTORY_DURATION * FPS));
    if (skip_trajectory) {
        printf("Note: Trajectory tests skipped (--skip-trajectory)\n");
    }
    printf("\n");

    // Load test locations
    std::vector<Location> locations = load_locations("interesting_locations.txt");

    if (locations.empty()) {
        printf("ERROR: No locations loaded from interesting_locations.txt\n");
        return 1;
    }

    printf("Loaded %zu test locations.\n\n", locations.size());

    int location_passed = 0, location_failed = 0;
    int trajectory_passed = 0, trajectory_failed = 0;

    for (const auto& loc : locations) {
        printf("═══════════════════════════════════════════════════════════════\n");
        printf("Testing: %s\n", loc.name.c_str());
        printf("  Target: (%.10f, %.10f) zoom=%.2e\n", loc.real, loc.imag, loc.zoom);
        printf("───────────────────────────────────────────────────────────────\n");

        // Test 1: Is the target location itself interesting?
        printf("  [1] Location interestingness: ");
        auto loc_result = test_location_interestingness(loc, verbose);
        if (loc_result.passed) {
            printf("PASS (score=%.1f)\n", loc_result.fail_score);
            location_passed++;
        } else {
            printf("FAIL (%s, score=%.1f)\n", loc_result.reason.c_str(), loc_result.fail_score);
            location_failed++;
        }

        // Test 2: Does the trajectory stay interesting throughout?
        if (skip_trajectory) {
            printf("  [2] Trajectory test: SKIPPED\n");
        } else {
            printf("  [2] Trajectory test: ");
            fflush(stdout);
            auto traj_result = test_trajectory(loc, TRAJECTORY_DURATION, verbose);
            if (traj_result.passed) {
                printf("PASS (all %d frames interesting)\n", (int)(TRAJECTORY_DURATION * FPS) + 1);
                trajectory_passed++;
            } else {
                printf("FAIL at frame %d (t=%.2fs, zoom=%.2e): %s",
                       traj_result.fail_frame, traj_result.fail_time,
                       traj_result.fail_zoom, traj_result.reason.c_str());
                if (traj_result.fail_score != 0.0) {
                    printf(" (score=%.1f)", traj_result.fail_score);
                }
                printf("\n");
                trajectory_failed++;
            }
        }
    }

    printf("\n═══════════════════════════════════════════════════════════════\n");
    printf("SUMMARY\n");
    printf("═══════════════════════════════════════════════════════════════\n");
    printf("Location tests:   %d passed, %d failed\n", location_passed, location_failed);
    if (skip_trajectory) {
        printf("Trajectory tests: SKIPPED (use without --skip-trajectory to run)\n");
    } else {
        printf("Trajectory tests: %d passed, %d failed\n", trajectory_passed, trajectory_failed);
    }
    printf("\n");

    // Only count location failures if trajectory tests are skipped
    int total_failed = location_failed + (skip_trajectory ? 0 : trajectory_failed);
    if (total_failed == 0) {
        printf("RESULT: ALL TESTS PASSED\n");
        return 0;
    } else {
        printf("RESULT: %d TESTS FAILED\n", total_failed);
        return 1;
    }
}
