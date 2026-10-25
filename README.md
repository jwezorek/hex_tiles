# Hex Tiles

A C++23 command-line tool for exploring substitution tilings built from five
colored tile types: yellow, blue, green, purple, and red. It generates finite
substitution patches and SVG diagrams of the substitution rules.

The broader goal is to investigate local matching rules and recognizability for
a potentially aperiodic hexagonal tiling system. The current implementation
generates and renders patches; it does not derive matching rules, prove
aperiodicity, or reconstruct full hexagons.

## Requirements

- A C++23-capable compiler
- CMake 3.16 or newer
- Eigen3 3.3 or newer
- Boost with a CMake package configuration
- Optional: Python 3 for SVG XML validation tests (standard library only)

The exact `triangle_grid` library has no Eigen or Boost dependencies, but
configuring the full project still requires both packages.

## Build

Run these commands from the project root, with the dependencies available to
CMake:

```sh
cmake -S . -B out/build -DCMAKE_BUILD_TYPE=Release
cmake --build out/build --config Release
```

If CMake cannot locate the dependencies, supply their installation prefixes
through `CMAKE_PREFIX_PATH`, or set `Eigen3_DIR` and `Boost_DIR` to their CMake
package directories.

For single-configuration generators such as Ninja, the executable is normally
`out/build/hex_tiles` (`hex_tiles.exe` on Windows). Multi-configuration generators
such as Visual Studio normally place it in `out/build/Release/`.
`CMakeSettings.json` also supplies Visual Studio x64 Debug and Release
configurations using Ninja.

## Usage

These examples assume the executable's directory is on `PATH`; otherwise, use
the path to the executable produced by your build.

```sh
# Default: geometric mode, golden lengths, green seed, depth 4
hex_tiles

# Geometric patch with regular leaf dimensions
hex_tiles --preset regular --seed green --depth 4 --out-dir output/regular

# Exact triangle-grid patch, including direction arrows
hex_tiles --mode combinatorial --preset regular --seed green --depth 4 --arrows --out-dir output/grid

# Geometric patch with custom leaf dimensions
hex_tiles --long 3 --short 1 --seed purple --depth 3 --out-dir output/custom

hex_tiles --help
```

Each run creates the output directory if necessary and writes:

- `patch.svg`: the generated patch for the selected seed and depth.
- `rules.svg`: the one-step substitution diagrams for all five tile types.

Open these files in a web browser or SVG viewer. Reusing an output directory
overwrites its generated SVGs.

### Options

| Option | Values / behavior | Default |
| --- | --- | --- |
| `--mode` | `geometric` or `combinatorial` | `geometric` |
| `--preset` | `golden` or `regular` | `golden` |
| `--long L --short S` | Custom finite leaf dimensions with `L > S > 0`; geometric mode only | Unset |
| `--seed` | `yellow`, `blue`, `green`, `purple`, or `red` | `green` |
| `--depth` | Nonnegative integer; 0 produces one canonical leaf | `4` |
| `--out-dir` | Directory for both SVG files | `output` |
| `--arrows` | Show direction arrows | Off |
| `--no-arrows` | Disable direction arrows | Off |
| `--help` | Print usage information and exit | — |

Supply custom lengths together, and do not combine them with `--preset`.
Combinatorial mode requires an explicit `--preset regular`; it rejects golden
and custom lengths. Invalid arguments and generation or output failures return
a nonzero exit code.

## How it works

**Geometric mode** uses Eigen transformations and parameterized leaf shapes.
The golden preset uses `(L, S) = (phi, 1/phi)`, where
`phi = (1 + sqrt(5)) / 2`; the regular preset uses `(2, 1)`.
At each substitution level, dimensions follow:

```text
L_next = 3 * L - S
S_next = L
```

Depth counts substitutions while keeping the chosen leaf dimensions fixed.
Higher levels use the corresponding level-specific placement rules; arbitrary
custom dimensions do not have a universal inflation scale.

**Combinatorial mode** represents regular patches with exact integer coordinates
on a triangular lattice. It records cell ownership, adjacency, and the
substitution hierarchy, including parent IDs and child slots. Euclidean
projection is used for SVG rendering, while the core uses checked integer
arithmetic.

Patch sizes grow rapidly with depth. The geometric generator caps total cached
leaf instances at 1,000,000; the combinatorial generator caps elementary cells
at 1,000,000. Oversized requests are rejected.

## Tests

CTest is enabled by default. To configure explicitly, build, and run the suite:

```sh
cmake -S . -B out/build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build out/build --config Release
ctest --test-dir out/build -C Release --output-on-failure
```

The tests cover lattice transformations, adjacency and overflow, substitution
coverage and hierarchy, comparison with geometric patches, and CLI validation.
If Python 3 is found during configuration, CTest also validates the generated
SVG XML and compares rendering geometry. Test output stays in the build
directory. Use `-DPython3_EXECUTABLE=<path>` when Python is not on `PATH`.

## Project layout

| Path | Contents |
| --- | --- |
| `src/main.cpp` | CLI parsing and SVG generation entry point |
| `src/hex.*`, `src/tile.*`, `src/geometry.*` | Geometric shapes, substitutions, and rendering |
| `src/triangle_grid.*` | Exact lattice primitives and transformations |
| `src/combinatorial_hex.*` | Exact substitutions, hierarchy, and cell ownership |
| `src/combinatorial_svg.*` | SVG rendering for exact patches |
| `tests/` | C++ tests, CLI checks, and Python SVG validation |
| `docs/` | Detailed design and implementation notes |

Further reading:

- [Parameterized substitution](docs/parameterized-substitution.md)
- [Exact triangle-grid mode](docs/combinatorial-grid.md)
- [Research goals and matching-rule exploration](hex-tiles-project.txt)

The design notes include historical validation reports; use the commands above
to check your current checkout and environment.

## License

[MIT](LICENSE), copyright © 2026 joe wezorek.
