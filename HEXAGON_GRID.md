# Derived, oriented hexagon grid

The existing exact triangle-grid substitution is unchanged. The new
`jw::hexgrid::convert(const jw::grid::patch&)` builds hexagons **after**
substitution, using only integer triangle cells and the stored integer leaf
affine frames. It never compares or snaps Euclidean coordinates.

## Definition and lattice

A hexagon is the six elementary triangular cells incident at one triangle-grid
vertex. This is established by the unit (`L=2, S=1`) canonical regions: a
triangle-colored hexagon is six one-cell subtiles, while a blue, green, or
purple hexagon is two three-cell trapezoidal subtiles. All six cells must have
the same source color and a compatible directed/reflected leaf-arrow motif.

Let the original triangle-grid basis be `e_q=(1,0)`, `e_r=(0,1)` and its
projected vectors `(1,0)`, `(1/2,sqrt(3)/2)`. Hexagon centers occupy one of
three phases `triangle.q - triangle.r (mod 3)`. The converter chooses the phase
with the largest count of **complete, valid, decorated** hexagons, breaking ties
by the fewest inconsistent complete candidates. A positive exact tie is an
error (the phase is ambiguous). Patches too small to contain any complete
hexagon return an empty hex patch.

`hexgrid::coordinate` is **not** a triangle-grid `vertex`. The hex-grid axial
basis, expressed in triangle-grid coordinates, is

    U = (1,1),   V = (-1,2).

With the selected phase `p`, the triangle-grid vertex corresponding to hex
coordinate `(q,r)` is

    center(q,r) = (p,0) + q*U + r*V.

The inverse uses exact divisibility by three, and throws for an off-grid
vertex. Neighbor offsets, in cyclic order, are `(1,0)`, `(0,1)`, `(-1,1)`,
`(-1,0)`, `(0,-1)`, `(1,-1)`.

## Orientation and decoration

Every cell retains its source leaf ID and its owner's canonical `+q` arrow
frame. Its arrow direction is the image of `+q` under the integer affine map,
encoded `0..5` for counterclockwise 60-degree turns from right. The determinant
supplies the subtile reflection bit.

A hexagon's **six-sector decoration** is read counterclockwise from -150 degrees
(sectors centered at -150, -90, -30, 30, 90, 150 degrees). The canonical
signatures, as `(subtile arrow direction, reflected)` pairs, are:

- **Yellow/red**: `(4,1), (4,0), (0,1), (0,0), (2,1), (2,0)`.
- **Blue/green/purple**: `(0,1), (0,1), (0,1), (0,0), (0,0), (0,0)`.

The stored hex orientation is the smallest counterclockwise rotation `0..5`
which maps that type's canonical six-sector decoration to the observed one.
The SVG arrow points along that stored rotation (0=right, 1=up-right, etc.).
It is not inferred from coordinates or from colors alone.

The first motif has **120-degree rotational symmetry**: yellow and red only
have two distinguishable orientation classes (0 and 1), so rotating either by
120 degrees does not change its stored orientation. The trapezoid motif has
six distinguishable rotations. Reflections of **both** canonical motifs are
rotation-equivalent (the trapezoid motif reflects to its 180-degree rotation;
the triangle motif reflects to its odd-rotation class). Thus no independent
hexagon-level reflection bit is needed for these five prototiles; the **leaf**
reflection states remain recorded through the source IDs. This conclusion
would need revisiting if additional asymmetric hex motifs are introduced.

## Completeness and error handling

For a selected phase, every elementary triangle has **exactly one** center
vertex on that phase. The center's six incident triangles consequently form
a partition; no triangle can belong to two hexagons. The converter verifies
this and stores the six source cells and IDs in each accepted `hex_tile`.

Six occupied cells with mixed types, incorrect subtile ownership, or a motif
which does not match the canonical six-sector decoration are **errors**.
Incomplete hexes on the edge of the generated patch are omitted and counted
in `hex_patch::rejected_boundary`. A missing cell strictly inside the root's
expected exact integer support (including all three neighboring cells within
the support) is an error. This distinguishes missing interior cells from normal
partial hexes along the outside. Removing cells directly from the patch edge
is treated as boundary clipping. The root support is derived from the root
node's type, level, and integer affine transform.

`hex_patch` has a sorted map keyed by `coordinate`, `insert`, `find`, `at`,
`size`, `bounds`, `placement.center/locate`, and the free `neighbors` operation.
No substitution rules are added for hexagons.

## Generate a hexagon SVG

With the existing CLI, use:

```sh
hex_tiles --mode combinatorial --preset regular --seed green --depth 4 --hex-svg --out-dir output
```

This generates the existing `patch.svg` and `rules.svg` **plus**
`hexagons.svg`. Every complete hexagon has a colored six-sided outline and a
black center arrow. Each SVG group also carries `data-type`, `data-q`,
`data-r`, and `data-orientation`. `--arrows` and `--no-arrows` continue to
control the old subtile visualization, not the mandatory hexagon arrows.
For example, the green depth-4 patch converts to **1,089 complete hexagons**
and rejects **77 boundary fragments**. The SVG output for a depth with no
complete hexagons is a valid, empty SVG.

## Tests

`tests/hex_grid_tests.cpp` uses only the Eigen-free grid library. It checks
coordinate arithmetic, exact inverse mapping and rejection, known counts for
all five seeds at depths 0–4, disjoint source-cell membership, duplicate
insertion, bounds, omission of a clipped exterior hex, rejection of an
interior hole and a mixed hex, preservation through all 12 dihedral transforms,
and SVG markers. The existing CLI test also checks `--hex-svg` output, and its
Python XML validator checks arrow directions and polygon coordinates.

The full executable requires the project's existing Eigen, Boost and C++23
configuration. The hex-grid core itself is independent of Eigen and Boost.
