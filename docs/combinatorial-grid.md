# Exact triangle-grid mode

The new mode implements the user's combinatorial substitution specification in
place, without changing the geometric implementation or committing changes.

## Design and alternatives

Vertices use signed 64-bit coordinates `(q,r)` in the basis `(1,0)`,
`(1/2,sqrt(3)/2)`. A cell is a rhombus anchor and an up/down bit. Up cells have
vertices `(q,r),(q+1,r),(q,r+1)`; down cells have vertices
`(q+1,r),(q+1,r+1),(q,r+1)`. Sorted endpoint pairs identify unoriented edges.
This is more compact than storing three vertices per cell; cube coordinates
would add a redundant coordinate. All arithmetic in the core is checked integer
arithmetic. Euclidean projection belongs only to the renderer and comparison tests.

The positive 60-degree rotation is `[[0,-1],[1,1]]`; reflection across the
vertical Euclidean axis is `[[-1,-1],[0,1]]`. Affine composition applies the
right operand first. Rule anchors are translated symbolically from the existing
formulas, using `(x,y)=(q+r/2, sqrt(3)*r/2)`. Reflected rule anchors are already
in the reflected frame, just as in the geometric rules.

Rules at level k use lengths obtained after k recurrences, starting at (2,1).
A top-down traversal uses rules at level depth-1 and composes parent * child.
Preorder occurrence IDs, parent IDs, child slots, levels and transforms retain
the hierarchy. Only leaves own cells. A cell-to-owner map is authoritative and
rejects duplicate occupancy. A generated edge incidence map supports adjacency,
subtile seams and patch perimeter queries without tolerances.

SVG outlines are traced from uncancelled directed cell edges, avoiding seams
inside a subtile. Discrete boundary edge labels preserve the geometric renderer's
full-tile edge emphasis (which differs from the patch perimeter). Arrows use
transformed lattice endpoints, with the reflected normal handled at rendering.
Rules diagrams use the same five-row layout as the geometric output.

## Implementation plan and validation

1. Add lattice primitives and tests for the twelve lattice symmetries, inverses,
   adjacency and arithmetic overflow.
2. Add independent ordered rules, exact row enumeration of canonical regions,
   hierarchy construction and ownership. Test child sequences, depth-one counts
   9/21/21/16/4, coverage and recursive leaf counts at depths 0,1,2,4.
3. Add grid-derived SVG rendering and CLI selection. Reject golden and custom
   lengths in combinatorial mode; preserve the geometric default.
4. Compare rules and ordered leaves against the geometric implementation only in
   tests. Build and run CTest, inspect SVG XML, and review the final changes.

Review focus: reflected arrow offsets; child-slot ancestry; oversized depths;
overflow at integer limits; invalid CLI modes and incompatible lengths.

## Usage

```text
hex_tiles --mode combinatorial --preset regular --seed green --depth 4 --out-dir output_grid
hex_tiles --mode geometric --preset regular --seed green --depth 4 --out-dir output_geometry
```

Both write `patch.svg` and `rules.svg`. Add `--arrows` to show directions;
`--no-arrows` disables them. The default remains geometric with golden lengths.
Combinatorial mode requires an explicit regular preset and rejects custom lengths,
including custom values numerically equal to (2,1).

`triangle_grid` is a separate CMake library with no geometric dependencies.
Its public API is in `triangle_grid.hpp`, `combinatorial_hex.hpp`, and
`combinatorial_svg.hpp`. The application and comparison tests still use the
project's existing Eigen/Boost setup.

```cpp
auto p = jw::grid::make_supertile(jw::grid::tile_type::green, 4);
for (const auto& [cell, owner] : p.occupancy) {
    const auto& leaf = p.nodes.at(owner);
    // leaf.parent and leaf.child_slot identify its immediate substitution role.
    // Follow parent IDs until no_parent to recover the complete ancestry.
    for (auto neighbor : cell.neighbors()) {
        auto found = p.occupancy.find(neighbor); // exact adjacency lookup
        // found == end: exterior; found->second != owner: subtile seam
    }
}
auto edges = p.edges(); // one incidence: perimeter; two: shared edge
```

IDs are deterministic within a generated patch, and remain valid as vector
indices. They are not persistent IDs across independently generated patches or
different depths. Patches are intended to be consumed as constructed; if edited,
ownership, leaf cell lists and hierarchy must be kept consistent.

## Validation and limits

Configure with `-DBUILD_TESTING=ON` (the default), build, then run:

```text
ctest --test-dir <build-directory> --output-on-failure
```

The C++ tests cover lattice symmetries, inverses, adjacency, overflow, ordered
rules, exact coverage at four dimension levels, recursion at depths 0/1/2/4,
ownership, ancestry, and comparison with the geometric implementation. The CLI
test covers every seed at those depths, incompatible parameters, large-depth
rejection, arrows, and preservation of the geometric default. If Python 3 is
available at configure time, an additional standard-library test parses the SVGs
and compares outlines, emphasized edges and arrow endpoints with the reference
SVGs. Set `-DPython3_EXECUTABLE=<path>` if Python is not on PATH. No Python package
installation is needed. CTest artifacts stay in the build directory.

Verified locally with MSVC 19.51 in Debug: all four CTest suites passed, including
90 parsed SVG files. The green depth-4 geometric and grid renders were also
visually inspected. Independent source review found no correctness issues.

Patches are limited to 1,000,000 elementary cells; oversized depths are rejected
before allocation/recursion. Occupancy uses an ordered map for simple exact
lookup, rather than a packed representation optimized for very large patches.
Arithmetic overflow throws instead of wrapping. Matching-rule discovery and
recognizability analysis are intentionally not implemented.

No Git metadata is present in this workspace. Changes were made in place;
no existing geometric source files were removed or rewritten, and nothing was
committed or pushed.
