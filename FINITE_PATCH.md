# Finite-patch substitution (exact combinatorial hexagons)

## Public API

Include `finite_patch.hpp` for:

```cpp
jw::grid::patch jw::hexgrid::to_triangles(const hex_patch&);
jw::grid::patch jw::hexgrid::substitute_triangles(const grid::patch&, unsigned steps);
jw::hexgrid::hex_patch jw::hexgrid::substitute(const hex_patch&, unsigned steps);
```

`to_triangles` and `convert` are declared in `hex_grid.hpp`. The latter now
accepts an optional explicit `lattice` to fix the hex phase:

```cpp
hex_patch convert(const grid::patch&, std::optional<lattice> placement = std::nullopt);
```

The triangles returned by `to_triangles` are a **leaf-only** `grid::patch`:
`nodes` and `leaves` each enumerate full base-length prototile occurrences;
`occupancy` records every elementary cell and its owning leaf. An original
supertile, parent ID, or leaf ID is not needed. Original IDs in `hex_tile` are
strictly diagnostic metadata and ignored when reconstructing the patch.

## Reconstructing original prototile occurrences

A hex's center follows the supplied exact `lattice` and axial coordinates. Its
six incident elementary cells are generated at that center. The existing
six-sector **motifs** provide the type, reflection and arrow direction of each
leaf. For yellow/red there are six unit triangular leaves; the other colors
have two three-cell trapezoidal leaves, in opposite three-sector groups.
Rotate the canonical sector motif by the hex's orientation (two distinct
orientations for yellow/red, six for the others). Recover each full dihedral
integer affine map by matching its canonical elementary-cell set to its
sector group, trying exact integer translations. A unique match is required.
All six sectors must be owned precisely once. No floating-point coordinates
or ancestry information are consulted.

## Placement under substitution

The source leaf prototiles have **base lengths** `(L,S)=(2,1)` regardless of
how many times the represented patch has already been substituted. A single
substitution replaces each such occurrence with the **existing**
`grid::rules({2,1})` children. These children still have base dimensions;
their *parent* has dimensions `(5,2) = next_lengths({2,1})`.

Crucially, **the parent translations cannot simply be copied from the old
occurrences**: lengths 1 and 2 at the regular-hexagon realization change to
different lengths depending on the subtile type and the original corner role.
For every original prototile polygon corner, identify its corresponding
corner of the new `(5,2)` parent polygon. For every vertex shared by two or
more original prototile corners, impose:

```
new_translation[b] - new_translation[a]
    = transformed_new_local_corner[a] - transformed_new_local_corner[b]
```

Here the linear dihedral part of each affine transform is unchanged. Solve
these exact integer potential-difference equations using the graph of
shared-corner incidences. The first leaf fixes the otherwise free global
translation. Reject inconsistent equations, overlaps of expanded parents,
changes to their shared elementary-edge adjacency pairs, and changes in
Euler characteristic of the occupied triangle complex (which would expose
unexpected holes or topology changes). Expand the
parent polygons using the existing child rules and verify that the child
union covers them without overlapping. The global mesh is thereby deformed
*coherently*, not uniformly scaled and not independently translated.

Repeat on the **new base-length children** for subsequent steps, resolving
the shared-corner graph again at every step. This use of the same base-level
rule for each generation is valid only because each intermediate patch is
represented at base leaf dimensions. It is **not** repeated application of
fixed rule transforms to the original locations. In particular the global
parent dimensions evolve via `(L,S) -> (3L-S,L)`; no constant similarity
ratio is assumed. Exact comparisons with deeper supertiles test this
behavior through two steps.

## Hex output and finite boundaries

The output is recognized using the original six-sector motif recognizer.
To fix its global hex-grid phase without voting across a possibly tiny
finite patch, independently substitute the first input hex (in map order)
through the requested depth. It shares the same first-leaf translation gauge
as the complete input. Try all three *exact* phase classes on that small
witness and require **exactly one** valid phase admitting a complete motif.
Apply that phase explicitly to the full output; otherwise report ambiguity.

The converter groups every occupied triangle according to this lattice.
Incomplete six-sector rings on the finite patch boundary are dropped and
counted in `rejected_boundary`. Complete but invalid decorated rings are
errors. Finite leaf-only triangle patches have their *entire* intended
support verified against their original occurrences before conversion:
missing or overlapping elementary cells cannot be silently discarded.
Fully absent hexes are allowed as intentional holes in an arbitrary finite
hex input; they form additional boundaries. A rooted supertile additionally
provides its enclosing polygon for detecting missing interior cells.

## Admissible inputs and limitations

A hex patch can be a connected nonconvex set, including patches with holes.
Every hex has a valid type, canonical orientation, and axial coordinate;
other `hex_tile` metadata may be absent. Reconstruction of all five hex
motifs is unique, but **not every arbitrary juxtaposition of colored/oriented
hexes is substitutable**: corner equations and shared-edge constraints must
be geometrically consistent. The method detects and rejects violations.

Disconnected triangular components are currently rejected for positive
substitution depth. Their relative translations after changing the length
parameters are not determined by the supplied tile labels and orientations;
arbitrarily retaining their old separation could conceal an overlap or
invent a registration. Empty patches and zero-step substitutions are valid.
All operations use checked triangle-grid arithmetic and retain the existing
1,000,000-cell limit. The current algorithm builds maps of triangle cells and
shared-corner constraints; it is intended for finite local patches, not
unbounded substitutions.

## Command-line demonstration

```text
hex_tiles --mode combinatorial --preset regular --seed blue --depth 4 \
    --finite-patch --patch-radius 2 --substitute-steps 2 --hex-svg \
    --out-dir output
```

The executable generates a supertile and its hex representation, chooses an
interior **complete** radius-2 axial hex ball (19 hexes), and calls the
provenance-free `substitute` API on **only those 19 hexes**. The outputs are:

- `patch.svg`, `rules.svg`: existing behavior
- `hexagons.svg`: existing optional complete source hexes
- `finite_original.svg`: extracted radius-2 hex ball and its orientation arrows
- `finite_substituted.svg`: hexes surviving two substitution steps, with arrows

Use `--substitute-steps 0`, `1` or `2`, or another nonnegative number;
`--patch-radius` also supports `0` for one hex. Increase source `--depth`
if no complete ball of the requested radius fits inside the supertile.
No original CLI behavior is changed unless a finite-patch option is selected.

## Tests

`finite_patch_tests` exercises every legal orientation of all five colors,
zero/one/two steps, reconstruction of exact leaf frames from supertiles,
independent deeper-supertiling region comparisons after axial registration,
connected mixed patches, deliberate holes, deterministic insertion order,
misstated provenance, overlap/missing-cell rejection, partial boundary
rings, disconnected inputs, and empty patches. The existing triangle-grid
and hex-grid tests remain applicable.
