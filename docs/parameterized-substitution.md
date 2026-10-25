# Parameterized hex-tile substitution demonstrator

This is a geometric C++23/Eigen demonstrator. The five labels, colors, arrow
indices, external-edge metadata, and ordered substitution rules are unchanged.
It does **not** reconstruct full hexagons or derive matching rules.

## Dimensions and presets

`jw::hex_lengths{L,S}` specifies **leaf** dimensions:

- `L`: long trapezoid base
- `S`: historically named small length
- `M = L - S`: remaining/middle length

The required domain is finite `L > S > 0`; `S` is allowed to exceed `M`.
The **golden** preset is `(L,S) = (phi, 1/phi)` and the **regular** preset
is `(2,1)`. All five canonical leaf shapes are constructed with those exact
lengths. For custom input, `M` is **not** silently normalized to 1.

In one substitution, leaf parameters `(L,S)` produce a parent at:

```
L_next = 3 L - S
S_next = L
M_next = L_next - S_next = 2 L - S
```

The same recurrence in `(M,S)` coordinates is `(2 M + S, M + S)`.
For regular input, `(L,S,M)` at levels 0–4 is:

```
0:  2,  1,  1
1:  5,  2,  3
2: 13,  5,  8
3: 34, 13, 21
4: 89, 34, 55
```

There is no universal inflation scale. `get_hex_rules(lengths)` returns the
**raw rigid** child placements at these leaf dimensions. Its parent geometry
is at `next_lengths(lengths)`; its children are at `lengths`.

## Depth and matrix conventions

`make_supertile(seed_type, depth, leaf_lengths)` returns leaf instances,
with one identity-transformed seed for depth 0. Depth 1 is one raw rule, and
depth *n* tiles the canonical parent shape at the *n*-times-updated lengths.
The original leaf shape dictionary is used to render every depth.

For each level, the generator caches the expanded patch for all five types;
for every parent and child in source order it appends all leaf instances in
that child's previous-level patch, using the transform `child.trans * leaf.trans`.
This is bottom-up composition, not repeated `inflate()` calls with fixed rules.
The low-level `inflate` helper is still available for single steps, but must
not be applied repeatedly with one unchanging raw rule dictionary except in
an independently justified construction. The generator caps total cached
leaf instances at 1,000,000 to avoid pathological depths.

In the golden case, the resulting depth-*n* patch is globally `phi^(2n)`
times the old normalized patch; after accounting for this global factor, the
same child sequence and placements are intended. The original program's
one-step patch followed by four expansions corresponds to **depth 5** here.

## Build and run

The existing CMake project remains C++23, Eigen3, and the existing Boost
configuration. Configure it in an environment with these dependencies:

```sh
cmake -S . -B out/build -DCMAKE_BUILD_TYPE=Release
cmake --build out/build --config Release
```

Run the resulting `hex_tiles` executable (on Windows the executable may be
under `out/build/Release/hex_tiles.exe` depending on the generator):

```sh
hex_tiles --preset golden  --seed green --depth 4 --out-dir output/golden
hex_tiles --preset regular --seed green --depth 4 --out-dir output/regular
hex_tiles --long 2 --short 1 --seed green --depth 3 --out-dir output/custom
hex_tiles --preset regular --seed green --depth 4 --out-dir output/regular --arrows
hex_tiles --help
```

`golden` is the default preset, `green` the default seed, and 4 the default
depth. `--long` and `--short` must be supplied together and conflict with
an explicit `--preset` in either order. Valid seeds: `yellow`, `blue`,
`green`, `purple`, `red`. Depth is a nonnegative integer. Invalid values,
excessive patch sizes, nonfinite derived geometry, and output I/O failures
produce an error and a nonzero exit code.

`--out-dir` is created if necessary; it receives **patch.svg** and
**rules.svg**. The rules image places each canonical leaf shape beside its
**unscaled** one-step child patch, in consistent row order, using measured
bounds and a shared drawing scale. SVG pixel sizing/viewBox is display-only.
Both SVGs draw exterior-edge outlines; arrows are disabled by default and
can be enabled with `--arrows` (or disabled with `--no-arrows`). Arrow inset
distances are proportional to their own segment length.

## Verification and limitations

The implementation was audited against the supplied placement tables and
source call sites. A separate Python numerical transcription of the actual
`place(...)` rows checked one-step polygon containment, non-overlap, and
area for all five types for `(M,S)=(1,1/phi), (1,1), (1,0.3), (1,1.8)`.
It also verified the child counts `(5,11,11,8,2)`, regular one-step areas
`(9,21,21,16,4)` in unit triangles, and a green regular depth-4 leaf count
of 3,571. That patch's triangle area was 6,765 and its maximum vertex
triangular-lattice coordinate deviation was about `1.5e-14`.

Eigen headers were unavailable. The three modified translation units
(`main.cpp`, `hex.cpp`, and `tile.cpp`) passed a GCC 14 C++23 **syntax-only**
check using a minimal declaration-only placeholder for Eigen; this is not a
build against real Eigen. **No executable was built or run and no SVGs were
visually inspected.** Run a local CMake build against Eigen and inspect both
SVG outputs; these checks do not prove that the real project compiles or
that the SVG output looks correct.
`to_full_tiles()` remains intentionally unimplemented.
