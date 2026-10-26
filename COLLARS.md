# Exhaustive hexagon collars

`jw::hexgrid::generate_collars` returns an exact, deduplicated catalog;
`jw::hexgrid::collars_to_svg` renders an existing catalog without enumeration or
mutation. Include `collars.hpp` and `hex_grid_svg.hpp` and link `triangle_grid`.

```cpp
#include "collars.hpp"
#include "hex_grid_svg.hpp"

int main() {
    const auto catalog = jw::hexgrid::generate_collars(6);
    jw::hexgrid::collars_to_svg("collars.svg", catalog, {5, 12.0});
}
```

A runnable version is `examples/collars.cpp`, built as `collar_example`:

```text
cmake --build out/build/collars-release --target collar_example
out/build/collars-release/collar_example.exe 6 collars.svg
```

The executable takes an optional radius and output path. Use the equivalent
executable path for your build configuration or platform.

`collar_catalog` contains `radius`, `inflation_depth`, `depth_basis`, optional
`coverage_radius`, and the sorted `collars`. Each collar is a vector of
`collar_hex { at, type, orientation }`. Its root is axial `(0,0)` and its support
is exactly the ball with `1+3*n*(n+1)` positions. Semantic keys omit ownership
IDs, lattice phase, observed sectors, and generation order. Coordinates are
relative to the root; orientations are parity for yellow/red and sixth turns
for blue/green/purple. Indices start at zero.

`canonicalize_collar(input, center)` translates the given center to zero, rotates
coordinates by `(q,r) -> (-r,q+r)` and orientations by one sixth turn, then takes
the lexicographically least sorted vector over all six proper rotations. It
rejects empty inputs, duplicate positions, invalid labels, and arithmetic
overflow. It never identifies reflections. It can also canonicalize non-ball
finite patterns; `generate_collars` supplies complete balls.

## Completeness and depth policy

The compiled constants in `src/collar_certificates.cpp` are the exact 19 legal
vertex seeds and 29 rooted radius-2 contexts from `research/collar-contexts.txt`.
Their provenance is the finite upper/lower and bootstrap certificates in
`research/collar-depth.md`, sections 4–5. Tests compare every compiled entry
with that fixture and check the 18 induced edge types. Runtime never reads
research files or runs Stage-1 discovery.

For every vertex seed, enumeration reconstructs the actual base triangular
prototiles, applies the repository's coherent shared-corner substitution to the
whole triangular patch through all selected steps, and only then recovers all
supported hexagons. Every complete radius-n ball is extracted, canonicalized,
and inserted into a sorted semantic set. No intermediate hex conversion clips
the triangular information. The core uses checked integer lattice operations;
it has no Eigen, Boost, renderer, or floating-point dependency.

The default depths for radii 1–6 are the **certified minimum catalog depths**
`1,2,2,3,3,3`. Counts are `14,29,50,77,110,149`; radius zero gives five colors at
depth zero. The catalog certificates compare complete sets with sufficient
coverage depths `1,2,3,3,4,4`. A catalog certificate does **not** assert geometric
coverage: `coverage_radius` is absent and `depth_basis` is `certified_catalog`.

Larger radii use the report's exact finite coverage procedure. At each candidate
depth, solve the corner equations at the actual Fibonacci lengths, obtain the
output phase from one descendant path, and form each source face's supported
output footprint by exact elementary-sector centroid tests. For **every**
three-face star in each of the 29 contexts, including stars not incident to the
central face, compute an inward graph-distance transform. Maximize that depth
at each output center over all stars, then minimize over all centers in the
expanded central face and all contexts. This is R(k), not a largest ball found
in one star. The certificate reproduces `R(0..4)=0,1,2,4,11`.

Search depths in increasing order until `R(k)>=n`. The checked length recurrence
and `S_k>4*n` supply a finite mathematical search bound; resource caps below may
stop the search with an exception. There is no monotonicity assumption, sampled
supertile enumeration, stabilization criterion, or extrapolated count formula.
These results have `depth_basis=coverage_certificate` and a `coverage_radius`.
Set `collar_options::use_coverage_depth=true` to request this policy also for
radii 1–6. Successful results are complete; exceptions never return a truncated
catalog. Coverage values are cached in process after exact computation.

## Supported partial reconstruction

The existing `convert(triangles, optional_placement)` behavior remains the
default: discard incomplete boundary hexagons and preserve the original root
boundary conventions. The new explicit mode is:

```cpp
auto supported = jw::hexgrid::convert(
    triangles, std::nullopt, jw::hexgrid::reconstruction::supported_sectors);
```

It validates complete base prototile occurrences, exact owner frames,
occupancy, conflicting labels/decorations, and phase agreement. A base
prototile's canonical hex center is `(0,0)` for yellow/red and `(1,0)` for the
other colors. Those exact frames determine the phase, even without a complete
six-sector motif. At that phase, the shared production motif logic uniquely
decodes any represented sector. An explicit placement must agree with every
owner. Wholly absent hexes remain absent, including intentional finite holes.

`hex_tile::observed_sectors` records the actual observed-sector mask. Only
observed slots of `triangles` contain provenance; unobserved slots of
`source_leaves` are `grid::no_parent`. Conversion does not fabricate observations
for the rest of the hexagon. `to_triangles` continues to reconstruct a whole hex
from semantic data, deliberately ignoring provenance. Existing complete-only
conversion now records the full mask 63 on its output hexes.

## Limits and failures

The exact implementation accepts radii at most **32** and inflation depths at
most **5**, with at most **1,000,000 triangle cells per seed**, **4,000,000 lattice
scan candidates per parent polygon**, and **1,000,000 stored semantic hex
labels per catalog**. The renderer also caps its input at 1,000,000 hex labels.
These are resource ceilings, not a promise that every radius through 32 fits.
`collar_options` can lower the depth and catalog-storage budgets; increasing
those options cannot bypass hard caps. Seed cell estimates are computed from
exact rule counts before substitution; ball and scan sizes are checked before
allocation. Coordinates, affine transformations, cross products, distance
increments, and length recurrences use checked 64-bit arithmetic. Depth/radius
limits bound all small loop counters and indexing arithmetic.

Unsupported resource requests throw `std::length_error`. Malformed labels and
SVG options throw `std::invalid_argument`; inconsistent observations throw a
logic error; checked coordinate operations throw on overflow. Use an optimized
build for coverage and larger catalogs: these are exhaustive finite procedures.
No cancellation or partial-result API is supplied.

The SVG options are column count (1–1000) and hex circumradius (4–1000 SVG units).
The renderer validates full ball support, uses the existing palette and
triangle-lattice projection, marks every orientation, and outlines each central
hexagon. Every panel has a stable index and each hex has axial, color,
orientation and panel attributes. Shared marker IDs are unique. Bounds include
margins for outlines, arrows, labels and title. An empty input produces an
explicitly labeled empty contact sheet. Open and write failures throw
`std::runtime_error`.

## Verification

CTest includes `collar_partial_reconstruction`, `collar_catalog_certificates`,
`collar_svg_generation`, and (when Python is available) `collar_svg_xml`.
The catalog test covers all counts 0–6, exact ball supports, rooted translation
and proper-rotation invariance, reflection distinction, sorted deterministic
keys, source data promotion, induced edge count, set equality at catalog and
coverage depths, R(0..4), radius-7 coverage selection, and resource errors.
Partial tests cover every decorated label and every complete owner fragment,
actual provenance, default compatibility, malformed support/frames/occupancy,
phase mismatches, conflicting observations and intentional holes. XML tests
inspect every panel, polygon, arrow, central outline, ID, orientation and layout
bound; generation tests check repeat output and catalog immutability.

Representative outputs are under the build directory's `collar-test-output`.
The radius-1 and radius-6 sheets were rasterized and visually inspected,
including a detailed radius-6 crop. All new tests pass. The full suite retains
the pre-existing Windows-only `hex_grid_reconstruction` cleanup failure: its
input stream is still open when it deletes the temporary SVG. That failure was
reproduced before implementation and left unchanged.

A further public-API smoke test generated radius 12 at computed coverage depth
5, producing 509 collars. This is a direct enumeration result, not a formula or
a claim of minimum catalog depth for that radius. The separate code review
reported no actionable findings.

## Changed files

- New exact API/core: `src/collars.hpp`, `src/collars.cpp`,
  `src/collars_detail.hpp`, `src/collar_certificates.cpp`, `src/collar_coverage.cpp`.
- Partial conversion and provenance: `src/hex_grid.hpp`, `src/hex_grid.cpp`.
- Independent catalog rendering: `src/hex_grid_svg.hpp`, `src/hex_grid_svg.cpp`.
- Tests: `tests/collar_partial_tests.cpp`, `tests/collars_tests.cpp`,
  `tests/collars_svg_tests.cpp`, `tests/collars_svg_validation.py`.
- Build/example/docs: `CMakeLists.txt`, `examples/collars.cpp`, `README.md`,
  `COLLARS.md`, `docs/plans/collars.md`.

The supplied research artifacts are unchanged. All work remains uncommitted in
the original checkout and branch. Final Release build succeeds; CTest is 9/10
passing, with only the reproduced baseline cleanup failure described above.
