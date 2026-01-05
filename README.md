# Mandelbrot Explorer

Ultra-fast console-based Mandelbrot fractal explorer with deep zoom capabilities using perturbation theory and double-double arithmetic.

![Mandelbrot Animation](assets/mandelbrot.gif)

## Features

- **Deep Zoom**: Explore to zoom levels beyond 10^30 using perturbation theory
- **Series Approximation (SA)**: 4-term polynomial skips up to 100% of iterations at deep zoom
- **Bilinear Approximation (BLA)**: Additional iteration skipping where SA is less effective
- **Block-Level Rendering**: Tile-based evaluation fills uniform regions without per-pixel computation
- **Reference Orbit Caching**: Avoids recomputation when panning at deep zoom
- **Double-Double Precision**: ~31 decimal digits of precision for reference orbit computation
- **AVX2 SIMD**: Optional 4x parallel pixel computation on supported CPUs
- **Smooth Animation**: Trajectory mode with ease-in-out cubic easing
- **Multiple Color Schemes**: 9 built-in palettes with rotation support
- **Real-time Navigation**: Pan, zoom, and rotate interactively
- **iTerm2 Image Mode**: High-resolution rendering using iTerm2 inline images

## Performance

Measured speedups with SA+BLA optimizations vs baseline (no SA):

| Zoom Level | Speedup | SA Skip |
|------------|---------|---------|
| 1e11 | **25x** | 100% |
| 1e12 | **28x** | 100% |
| 1e13 | **31x** | 100% |
| 1e14 | **31x** | 100% |

Even in challenging regions with low SA effectiveness (~5% skip), block-level rendering provides 2-3x speedup.

## Building

```bash
# Standard build (portable)
make

# With AVX2 optimization (faster on supported CPUs)
make avx2

# Native optimizations for current CPU
make native
```

## Usage

```bash
# Interactive mode
./mandelbrot

# Start at specific position and zoom
./mandelbrot --pos -0.7+0.3i --zoom 1e6

# Automatic exploration mode
./mandelbrot --auto

# Animated trajectory (zoom from default to target over 60 seconds)
./mandelbrot --pos -0.7115114743-0.3078112463i --zoom 1.86e+11 --auto=60
```

### Command Line Options

| Option | Description |
|--------|-------------|
| `--pos <re+imi>` | Target position (e.g., `-0.5+0.3i`) |
| `--zoom <value>` | Target zoom level (e.g., `1e6`) |
| `--angle <degrees>` | Target view angle (e.g., `45`) |
| `--auto [N]` | Auto exploration, or trajectory over N seconds |
| `--image [WxH]` | iTerm2 image mode (e.g., `--image=800x600`) |
| `--benchmark` | Compute one frame and print timing (no interactive mode) |
| `--no-sa` | Disable Series Approximation |
| `--no-bla` | Disable Bilinear Approximation |
| `--no-block` | Disable block-level tile filling |
| `--no-cache` | Disable reference orbit caching |
| `--width N` | Set render width (for --dump-iter) |
| `--height N` | Set render height (for --dump-iter) |
| `--dump-iter PATH` | Dump iteration buffer to file (for testing) |
| `--help` | Show help message |

### Interactive Controls

| Key | Action |
|-----|--------|
| Arrow Keys | Pan view |
| SHIFT + Up/Down | Zoom in/out |
| SHIFT + Left/Right | Rotate view |
| Z (macOS) | Toggle arrow mode (Pan <-> Zoom/Rotate) |
| C/V | Rotate color palette |
| 1-9 | Switch color schemes |
| +/- | Adjust max iterations |
| I | Toggle iTerm2 image mode |
| R | Reset view |
| Q/ESC | Quit |

## Technical Details

### Series Approximation (SA)

At deep zoom, most pixels follow nearly identical iteration paths. Series Approximation exploits this by computing polynomial coefficients that approximate the perturbation:

```
δZ_n ≈ A_n * δC + B_n * δC² + C_n * δC³ + D_n * δC⁴
```

The coefficients follow recurrence relations derived from the Mandelbrot iteration:
- `A_{n+1} = 2*Z_n*A_n + 1`
- `B_{n+1} = 2*Z_n*B_n + A_n²`
- `C_{n+1} = 2*Z_n*C_n + 2*A_n*B_n`
- `D_{n+1} = 2*Z_n*D_n + 2*A_n*C_n + B_n²`

The polynomial is evaluated using Horner form for efficiency: `(((D*δC + C)*δC + B)*δC + A)*δC`

At zoom 1e12, SA typically skips 100% of iterations, providing massive speedup.

The validity check ensures the approximation error stays below tolerance:
```
|D_n|² * |δC|⁶ < ε² * |A_n|²
```

**Automatic threshold**: SA is automatically disabled at zoom > 1e12 where the 4-term polynomial approximation loses precision. BLA continues to provide iteration skipping at any zoom level.

### Bilinear Approximation (BLA)

BLA provides additional iteration skipping where SA is less effective (near minibrots). Unlike SA which approximates in terms of δC only, BLA is linear in both δz and δC:

```
δz_{m+l} = A_l * δz_m + B_l * δC
```

BLA coefficients are organized in a hierarchical table:
- Level 0: 1-iteration BLAs (A = 2Z, B = 1)
- Level 1: 2-iteration BLAs (merged from level 0)
- Level k: 2^k-iteration BLAs

Merging formula: `A_{merged} = A_y * A_x`, `B_{merged} = A_y * B_x + B_y`

### Block-Level Rendering

Instead of computing each pixel individually, the renderer first attempts to fill 16x16 tiles:
1. Sample the 4 corners of each tile
2. If all corners agree (all escaped or all bounded with similar iteration counts), fill the entire tile
3. Otherwise, fall back to per-pixel computation

This provides significant speedup in uniform regions (interior of cardioid, exterior far from boundary).

**Automatic threshold**: Block-level filling is automatically disabled at extreme zoom (pixel size < 1e-15) where corner interpolation can introduce visible artifacts.

### Reference Orbit Caching

The reference orbit is cached and reused when:
- Center position hasn't changed
- Max iterations haven't increased
- SA settings are the same

This allows smooth panning at deep zoom without recomputing the expensive reference orbit.

### Perturbation Theory

At deep zoom levels, standard double-precision floating point loses accuracy. This explorer uses perturbation theory:

1. Compute a high-precision reference orbit using double-double arithmetic
2. For each pixel, compute only the small perturbation from the reference
3. Detect and recover from "glitched" pixels where perturbation becomes inaccurate

### Double-Double Arithmetic

Double-double uses two `double` values to achieve ~31 decimal digits of precision (vs ~16 for standard double). This is sufficient for zoom levels up to approximately 10^30.

### AVX2 Optimization

When built with `-mavx2`, the explorer processes 4 pixels simultaneously using SIMD instructions, providing significant speedup on modern CPUs.

### iTerm2 Image Mode

When running in iTerm2 (detected via `LC_TERMINAL` or `ITERM_SESSION_ID`), the explorer can render using iTerm2's inline image protocol. This provides much higher resolution output (default 640x400 pixels) compared to terminal character cells. The image is encoded as PPM and transmitted via base64.

## Testing

The project includes comprehensive tests for optimization correctness:

```bash
# Run all tests
make test-all

# Run specific test suites
make test              # Perturbation vs direct DD comparison
make test-sa           # Series Approximation polynomial evaluation
make test-trajectory   # Animation frame interestingness
make test-flags        # Optimization flag combinations (64x40)

# Run flag tests at higher resolution
./test_flags --highres   # 256x160
./test_flags --fullres   # 640x400
```

The flag matrix tests verify that all optimization combinations produce correct results:
- Compares optimized rendering against fully-disabled baseline
- Tests 8 locations at various zoom levels (1e6 to 1e18)
- Validates escape/bounded classification and smooth iteration values
- Detects block-level tile artifacts

## Requirements

- C++17 compiler (clang++ or g++)
- POSIX terminal with ANSI escape code support
- Optional: CPU with AVX2 for SIMD acceleration
- Optional: iTerm2 for high-resolution image mode

## License

MIT License
