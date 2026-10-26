# Radius-1 periodic counterexample using all five colors

A periodic tiling exists satisfying every legal radius-1 collar and using every
color. Period vectors in axial coordinates are (7,0) and (0,7). The 49-position
cell contains yellow 12, blue 9, green 3, purple 24, red 1, with the orientations
recorded in out/build/radius1-periodic.txt and radius1-periodic.json.

The dependency-free research/radius1_periodic_search.py reads the generated
radius-1 SVG's semantic metadata, expands all rotations (76 distinct oriented
neighborhoods), and searches Hermite-normal-form tori using bitset constraint
propagation and exhaustive backtracking. It requires every color and fixes red
orientation zero at the origin by translation/proper-rotation symmetry. Every
reported solution is independently checked by direct wrapped neighborhood
membership rather than the search's overlap tables.

The independent research/radius1_periodic_verify.cpp reads the witness and
compares all 49 wrapped neighborhoods with canonicalize_collar against a fresh
generate_collars(1) production catalog. All pass. It also verifies all five
colors, and renders the cell and a 3x3 repetition using the existing renderer.
The witness uses collar indices 1,3,4,5,6,7,8,9,10,11,12,13 (zero based).

Reproduction from the repository root, with the bundled Python or any Python 3:

    python research/radius1_periodic_search.py --area 49 --seconds 1

Build the verifier from an MSVC developer shell after building triangle_grid:

    cl /nologo /std:c++latest /EHsc /O2 /MD /Isrc research/radius1_periodic_verify.cpp out/build/collars-release/triangle_grid.lib /Fe:out/build/radius1_periodic_verify.exe /Fo:out/build/
    out/build/radius1_periodic_verify.exe out/build/radius1-periodic.txt out/build/radius1-periodic

The exact search completed all period-cell areas 5 through 48 with no unresolved
cases, and found the witness during area 49 (1917 cases processed in total).
Areas below 5 cannot contain all five colors. The affirmative counterexample
itself needs no claim about minimality or solver performance.

The proposed blue/red-only construction fails the full radius-1 rules: among
the 76 oriented collars, the only two having exclusively blue/red labels are
red-centered. There is no allowed blue-centered collar using just those colors.
Thus blue neighbors of a red must bring in additional colors. The 49-cell
witness establishes the desired all-color counterexample despite that simpler
construction failing.

No production source or build configuration was changed for this investigation.
