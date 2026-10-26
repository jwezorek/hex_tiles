# Inflation depth for exhaustive collars

Investigation of the current branch, 2026-10-25. Production sources were left
unchanged. No commit or archive was made.

**Result:** vertex seeds suffice. For the supported, partially reconstructed
hexagons defined below, the exact minimum coverage depths for radii 1 through 6
are **1, 2, 3, 3, 4, 4**. The exact minimum catalog depths are
**1, 2, 2, 3, 3, 3**. In particular, a radius-6 catalog needs three inflations,
whereas covering every output position by a particular inflated source star
needs four.

There is an exact finite procedure for `k_cover(n)` for every n, using just
**29 certified source radius-2 contexts**. A sufficient depth always exists;
it grows logarithmically in n. A closed expression for the exact minimum at
arbitrary n has not been proved.

The numerical assertions here are finite exhaustive computations together with
the reductions proved below, rather than inferences from catalog stabilization.

## 1. Precise reconstruction and legality conventions

The intended tiling space is the substitution hull: infinite tilings obtained
as interior limits of the repository's supertiles, with their colored,
oriented prototile occurrences. Equivalently one may use complete substitution
hierarchies and their limits. Its finite patches have predecessors in such a
hierarchy. Uniqueness of a predecessor (recognizability) is unnecessary here.

Predecessor existence for the interior-limit definition follows by taking
predecessors of larger and larger supertiles, translating the preimage of a
fixed point into a bounded window, and taking a subsequence. There are finitely
many local decorated types; the deformations below are proper, with bounded
distortion at each fixed step. Windows with boundaries escaping to infinity
therefore have predecessor windows with the same property. This also gives
surjectivity of substitution on this hull, up to translation. This convention
excludes arbitrary infinite arrangements that merely satisfy some local tests.

For a triangular patch P, use the following precise interpretation:

> `Hex(P)` consists of every output hexagon with at least one elementary sector
> present in P, with its uniquely decoded color and orientation, on the correct
> output center lattice.

This is exactly the supported reconstruction promised by the question. In the
actual motifs, **one decorated sector suffices** once the lattice phase is known.
For yellow/red, the two orientation parities give distinct decorations at each
fixed sector. For blue/green/purple, the observed arrow direction selects the
orientation, and its reflection bit checks the sector group. The observed leaf
also supplies the color. These are direct checks of `triangle_motif` and
`trapezoid_motif` in `src/hex_grid.cpp`.

Additional information needed is the owner's color and exact dihedral frame,
the sector position, and consistent coordinate registration/phase. Colors or
undecorated occupied cells alone are insufficient. Ancestry is unnecessary.
In fact a full base prototile already fixes a hex center: its canonical center
is `(0,0)` for yellow/red and `(1,0)` for the other three. Its affine frame thus
fixes the global phase without a majority vote or a complete six-sector ring.

An output hexagon intersecting a substituted parent in **positive area** has
at least one whole elementary sector in that parent: parent boundaries follow
elementary triangle edges. Thus all such hexagons are reconstructible. Merely
touching a parent at a point or along an edge does not supply a sector.

The exact minima in this report concern this supported reconstruction. If
`Hex` is instead allowed to deduce wholly absent hexagons by arbitrary global
legality reasoning, that is a stronger, unspecified inference operation; the
sufficiency results remain valid but the necessity statements would have to
be reconsidered. The research probe does not implement that operation.

## 2. Geometry derived from the repository

The relevant implementations are `grid::rules`, `grid::canonical_region`,
`hexgrid::to_triangles`, and the shared-corner solver in
`hexgrid::substitute_triangles`. The five objects called triangular prototiles
include one-cell triangles and three-cell trapezoids at base lengths; these
must not be replaced by five abstract uniformly scaled triangles.

Put `M=L-S`. After k inflations,

\[
 S_k=F_{2k+1},\qquad M_k=F_{2k+2},\qquad L_k=F_{2k+3},
\]

where `F_1=F_2=1`. Equivalently `(M,S) -> (2M+S,M+S)`.

| k | L | S | M |
|---:|---:|---:|---:|
| 0 | 2 | 1 | 1 |
| 1 | 5 | 2 | 3 |
| 2 | 13 | 5 | 8 |
| 3 | 34 | 13 | 21 |
| 4 | 89 | 34 | 55 |

Reconstruct an input hex into its six triangles or two trapezoids, change their
corner roles to these lengths, and match shared corners, just as the finite
patch solver does. The union of the substituted parents is:

* Yellow: a regular hexagon of side M.
* Red: a regular hexagon of side S.
* Blue/green: a centrally symmetric equiangular hexagon with cyclic side
  lengths `(M,S,M,M,S,M)`, in its appropriate rotated frame.
* Purple: the same construction with `(S,M,S,S,M,S)`.

For the last two cases, write h for the trapezoid height and b for `L-h`.
Before rotation/translation, the union of the two expanded trapezoids has
vertices

\[
 (0,0),(h,-h),(h+b,-h),(h+b,0),(b,h),(0,h).
\]

This follows from the canonical trapezoid corners
`(0,0),(L,0),(L-h,h),(0,h)` and the opposite reflected half.
Its center is `(L/2,0)`, potentially half-integral. This auxiliary deformed
source center is not generally an output hex center.

An inflated source vertex patch is the union of three such hexagons, joined
along their original three interfaces. Its outer boundary is the image of
the original star boundary; shared internal edges cancel. Each original
boundary edge has length S or M after inflation. The corner equations fix
registration up to one common translation. These facts explain both the
overlap and the variation between vertex types. No Euclidean similarity of
the regular output grid is involved.

When comparing an independently substituted finite star with the full tiling,
register them by a corresponding source corner and its expanded image. The
library's first-leaf translation gauge can differ between patches. Equation
(4) places all candidate stars in one shared corner gauge, so its inclusions
compare actual relative coordinates and identities, rather than copies found
elsewhere in the tiling. Triangular information is retained through all k
steps; there is no intermediate hex conversion that clips the patch.

The ordered rules have translations linear in M and S. Composing the rules at
successive length levels gives these same outer parent polygons: the role of
each outer corner is preserved, and the level recurrence is `(3L-S,L)`.
Internal children are the repository's actual rules. The probe separately
compares the direct level-dependent construction with repeated finite-patch
substitution for all 19 vertex types through depth 3.

### A useful metric bound

On triangle coordinates use

\[
 N(q,r)=\max(|q|,|r|,|q+r|).
\]

There is an auxiliary piecewise affine homeomorphism `f_k` from the regular
source honeycomb to the deformed parent honeycomb. Map each source hex center
to the center of its expanded polygon and map its six fan triangles affinely.
On common source edges the maps agree. This is a convenient map of the parent
domains, not a claim that elementary cells substitute affinely.

It satisfies

\[
 S_k N(p-q)\ \leq\ N(f_k(p)-f_k(q))\ \leq\ M_k N(p-q). \tag{1}
\]

Here is a direct check of the distortion, which is also independent of color.
Up to lattice symmetries, the nonregular fan matrices are

\[
 A=\begin{pmatrix}(h+b)/2&(b-h)/2\\0&h\end{pmatrix},\qquad
 B=\begin{pmatrix}b&(b-h)/2\\0&h\end{pmatrix},
\]

and their reflected copies. Evaluate each matrix and its inverse on the six
vertices of the N-unit ball. Their operator norms are at most `max(h,b)` and
`1/min(h,b)`, respectively. Regular fans are scalar maps. Integrating along
straight segments gives (1) globally. The maps are homeomorphisms because
the positive convex parent hexagons share matching edges and have three
120-degree angles at every vertex. The induced complete, flat plane mesh has
no folds or gaps.

## 3. Why a sufficient finite depth exists

In the unit source honeycomb, every open N-ball of radius less than `1/2`
is contained in some vertex star. One elementary way to check this is that
two disjoint unit triangular-lattice edges have N-distance at least 1.
Consequently the ball can only meet a collection of honeycomb edges that
pairwise share endpoints. The honeycomb has no 3-cycles, so those edges share
a common vertex. All faces reached by the ball belong to that vertex's star.
If it meets no edge or only one edge, select a vertex of the containing face
or an endpoint of that edge.

Output axial coordinates have triangle basis
`U=(1,1), V=(-1,2)`, as in `hex_grid.hpp`. Each graph step has N-length 2.
All centers in a radius-n output collar therefore have N-distance at most
`2n` from its center. By (1), their preimages lie within `2n/S_k` of one source
point. Hence

\[
 \boxed{S_k>4n\quad\Longrightarrow\quad\text{coverage at depth }k.} \tag{2}
\]

Each of those output centers lies in the selected expanded star, so at least
one of its incident elementary sectors is available. The assumed decoder
recovers it. This proves existence uniformly over all legal tilings and all
output positions, without assuming an exhaustive V.

There is also a simple necessary bound. A unit source vertex star has
N-diameter 4. Its deformed domain has diameter at most `4M_k`; supported output
centers extend at most 1 beyond it on either end. A radius-n output ball has
diameter `4n` in N, using opposite multiples of U. Therefore coverage requires
`4n <= 4M_k+2`, or, for integers, `n <= M_k`. Thus

\[
 \min\{k:M_k\geq n\}\ \leq\ k_{\rm cover}(n)
 \ \leq\ \min\{k:S_k>4n\}. \tag{3}
\]

These are conservative bounds; the exact computation below is stronger.
Both sequences grow at rate `(3+sqrt(5))/2`, so the depth is
`Theta(log(n+1))`, not proportional to n.

## 4. One-step vertex coverage and completeness of Stage 1

There is a stronger special fact for *vertex configurations* than for collars.
Take an output honeycomb vertex w. The six elementary triangles incident at w
belong to its three output hexagons, two sectors per hexagon. The expanded
source faces containing w are either one face, the two faces of an edge, or
the three faces of a vertex. They are all incident to some source vertex v.
Every elementary triangle at w is in one of those faces, so `Sigma(S(v))`
contains all six. It determines the three output hex identities.

This argument uses the actual edge-aligned parent mesh and the supported
decoder. **Every output vertex type comes from one source vertex type in a
single step.** Separate edge seeds or a larger seed are unnecessary for this
operation. Every input edge adjacency is part of the vertex configuration at
either endpoint, so E is also obtained by taking pairs in V.

This alone does **not** make a forward-closure computation from an arbitrary
seed complete. There can be additional invariant or nonextendible types.
The following independent finite upper/lower certificate resolves that issue
for the present rules. Details are in `collar_probe.cpp` mode `admissible`;
its output is saved in `collar-certificate.txt`.

There are 22 individual decorated hex labels: two orientations each for yellow
and red, six each for the other colors. Fix the three positions
`(0,0),(1,0),(0,1)`. All vertex configurations, including both vertex classes,
are represented by the `22^3=10648` assignments and their proper rotations.
Reflections are never quotiented out.

1. Exhaust all 10648 assignments. Require consistent shared-corner expansion
   and compatible partial decorated output motifs. There are **39** surviving
   rotation/translation classes U. All legal vertex types must be in U.
   Testing through two and three steps still leaves 39; mere substitutability
   is visibly insufficient.
2. Compute `U intersect VertexTypes(Sigma(U))`. By the one-step argument and
   predecessor existence, every legal vertex type survives. There are **20**.
   Repeating this pruning stops at 20, including a spurious type
   `yellow(0)@(0,0), green(2)@(1,-1), green(4)@(1,0)`.
3. Enumerate every radius-1 assignment whose six vertex triples lie in those
   20 types. This is a finite cycle constraint on the six neighbors; there
   are **20** classes. All legal radius-1 neighborhoods are among them; no
   global extendibility is assumed.
4. Substitute these 20 seven-hex patches once. Collect output vertex types
   only at vertices lying in the expanded **central source face**, including
   its boundary. That face's source radius-1 patch includes all source faces
   incident to any such point. Hence it supplies the six sectors used in the
   preceding vertex argument. All output positions are represented by choosing
   a containing source face. The resulting upper bound is **19** types.
5. Independently, **all 19 occur** in depth-4 supertiles of the five repository
   seed types. These are genuine legal occurrences, not just passing local
   tests: each of the five seed types occurs strictly inside the green depth-4
   supertile with positive-determinant frame. The verified leaf IDs are
   blue 4, purple 5, yellow 6, green 41, red 335. In particular green contains
   an interior green. Nesting these occurrences and inflating gives unbounded
   clearance from the outer boundary, since S grows without bound. Every
   witnessed finite supertile patch therefore extends to an infinite
   substitution-hull tiling, up to proper rotation.

The upper and lower sets are equal, establishing **exactly 19 legal vertex
types**. The probe also checks forward closure of that exact set. Its pair
subconfigurations give **exactly 18 edge types**. The sample counts at depths
4, 5 and 6 are all 19, but that observation is not the completeness proof.

## 5. An exact, fixed-size reduction for coverage at every depth

Here is the important bound on the source context. Suppose the output center x
lies in expanded source face A, and an elementary sector of x meets another
expanded source face B. Pull back a point of that sector lying in B. In source
triangle coordinates, the displacement between source centers A and B has
N-norm at most

\[
 1+1/S_k+1\leq 5/2\qquad(k\geq1).
\]

For any source axial displacement c,

\[
 d_{\rm hex}(0,c)\leq\tfrac23 N(c_qU+c_rV).
\]

Thus the integer source graph distance between A and B is at most 1. Any star
that supplies x must contain one such B. The other two faces of that star are
neighbors of B, at source distance at most 2 from A. **Every eligible source
star, including stars not incident to A, is therefore entirely present in
A's source radius-2 neighborhood.** For k=0 this is immediate directly.

This reduction does not require a source neighborhood with radius depending
on n or k. The output collar may extend farther; its membership in each
eligible star image is tested against that image's actual finite boundary.

### Certifying the 29 radius-2 contexts without assuming local extendibility

The following bootstrap uses supersets until coverage justifies an exact set.

| Calculation | Number of classes | Interpretation |
|---|---:|---|
| Radius-1 assignments compatible with the exact 19 V | 18 | Legal superset |
| Join those radius-1 templates at all seven centers of a source radius-1 ball | 55 | Legal radius-2 superset |
| Check one-step coverage on all 55 | minimum R(1) = 1 | Exhaustive lower certificate |
| Radius-1 collars from once-inflated legal V | 14 | Consequently exact legal radius-1 catalog |
| Join those exact 14 templates as above | 37 | Legal radius-2 superset |
| Check two-step coverage on all 37 | minimum R(2) = 2 | Exhaustive lower certificate |
| Radius-2 collars from twice-inflated legal V | 29 | Consequently exact legal radius-2 catalog |

The joins enumerate all six proper rotations of every rooted template and
require identical labels on every overlap. They are exhaustive because any
actual radius-2 neighborhood restricts to one of these templates around each
of its seven interior centers. Extra locally compatible assignments are
retained until a coverage certificate removes the need to assume they are
legal. This is why neither `18=14` nor `37=29` was assumed.

All 29 final contexts are independently witnessed in depth-4 supertiles and
have the infinite extensions established above. Their full data, along with
the 19 vertex seeds, are in `collar-contexts.txt`.

### Exact definition and algorithm for R(k)

For each of those 29 rooted contexts c:

1. Solve the shared-corner equations directly at `(L_k,S_k)`. Call the expanded
   polygon of source hex a `D_k(c,a)`. Determine the output lattice phase from
   one descendant base prototile. Only one child path of length k is needed.
2. For every source vertex v whose three faces are in c, compute
   `H_k(c,v)`: output centers having at least one elementary sector in the
   union of the three polygons. This is `Hex(Sigma^k(S(v)))`.
3. For x in the output center lattice inside the expanded **central face**,
   compute

\[
 \rho_k(c,v,x)=d_{\rm hex}(x,\Lambda_k\setminus H_k(c,v))-1.
\]

   Set this to -1 when x is absent. An inward breadth-first distance transform
   computes it exactly: boundary hexes start at zero. Then set

\[
 \boxed{R(k)=\min_{c}\ \min_{x\in\Lambda_k\cap D_k(c,0)}
                      \ \max_v\rho_k(c,v,x).} \tag{4}
\]

The **maximum over all eligible stars is essential**. Eroding just one star,
or taking its largest contained collar, would answer a different question.

This is an exact finite procedure. The source-context bound proves that no
outside star can improve the maximum. Every position in every legal tiling
has a containing source face and thus a represented case. Conversely every
context is legal and its finite neighborhood extends to an infinite tiling,
so a minimizing case is an actual obstruction, not an inadmissible pattern.

Polygon membership uses integer cross products of three-times triangle
centroids. Boundaries lie along elementary grid edges, so a centroid test
classifies the entire elementary triangle exactly. There are no floating
point tests. Colors/orientations enter the polygon shapes and phase; at a
legal overlap the reconstructed identities agree with the infinite tiling.

Finally,

\[
 \boxed{k_{\rm cover}(n)=\min\{k\geq0:R(k)\geq n\}.} \tag{5}
\]

Equations (4)-(5), the 29 certified contexts, and the corner rules give the
requested exact function by algorithm. For any given n, (2) bounds the search
to finitely many k. There is no need to presume monotonicity of R or use binary
search. A mathematical implementation can use arbitrary-size integers; the
present research executable uses the repository's checked 64-bit coordinates
and retains its resource limits. Its tested scope is through k=4.

## 6. Computed exact values and explicit obstructions

| k | R(k), exact | Range of minimum coverage radii among the 29 rooted contexts |
|---:|---:|---:|
| 0 | 0 | 0–0 |
| 1 | 1 | 1–1 |
| 2 | 2 | 2–3 |
| 3 | 4 | 4–7 |
| 4 | 11 | 11–18 |

Every entry is obtained by the exhaustive reduction above. A failure at
`R(k)+1` is recorded with its actual source patch, output phase, target, and
one missing output coordinate for **each** eligible source star in
`collar-geometry.txt`. Other stars fail to supply even the target.

For example the same legal source radius-2 context supplies obstructions at
all these depths. Its central hex is purple, orientation 4. Here is its full
source data, with columns `(q,r,type,orientation)`, using
`0=yellow, 1=blue, 2=green, 3=purple, 4=red`:

```text
-2, 0,3,0    -2, 1,0,0    -2, 2,1,2
-1,-1,2,4    -1, 0,2,2    -1, 1,0,0    -1, 2,3,0
 0,-2,3,2     0,-1,2,0     0, 0,3,4     0, 1,3,2     0, 2,3,4
 1,-2,0,0     1,-1,0,0     1, 0,3,0     1, 1,0,0
 2,-2,0,0     2,-1,2,4     2, 0,2,2
```

The direct corner gauge anchors the first source prototile origin, as in the
probe. Failure targets are:

| k | Output lattice origin, triangle coordinates | Output target, axial coordinates | Failing collar radius |
|---:|---|---|---:|
| 0 | (0,0) | (0,0) | 1 |
| 1 | (2,0) | (2,0) | 2 |
| 2 | (1,0) | (13,-2) | 3 |
| 3 | (2,0) | (38,-6) | 5 |
| 4 | (0,0) | (106,-18) | 12 |

At k=3 the six eligible source vertices, in **source triangle coordinates**,
and their missing output centers, in **output axial coordinates**, are:

| Source vertex | Missing hex center |
|---|---|
| (-1,0) | (40,-11) |
| (-1,1) | (33,-6) |
| (0,-1) | (35,-1) |
| (0,1) | (33,-6) |
| (1,-1) | (33,-6) |
| (1,0) | (33,-6) |

Each missing center is at graph distance 5 from `(38,-6)`. All source stars
outside this eligible set lack the target itself. This gives an explicit
legal radius-5 failure at depth 3, and hence also a radius-6 failure. At depth
0 a star has only three hexes, so it cannot contain a seven-hex radius-1 ball.
The analogous complete depth-1 and depth-2 witnesses are in the saved output.
Together these obstructions prove minimality, not just sufficiency.

## 7. Catalog depth is smaller

The research experiment evaluates the assumed partial reconstruction on the
19 legal vertex seeds, collects **every** complete radius-n collar, and
canonicalizes translations and the six proper rotations. It does not merge
reflections. The successive counts are:

| Inflation depth | n=1 | n=2 | n=3 | n=4 | n=5 | n=6 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| 1 | 14 | 7 | 0 | 0 | 0 | 0 |
| 2 | 14 | 29 | 50 | 49 | 16 | 0 |
| 3 | 14 | 29 | 50 | 77 | 110 | 149 |
| 4 | 14 | 29 | 50 | 77 | 110 | 149 |

Every extracted collar is legal because each seed has a legal infinite
extension. Depth 4 covers every position for all six radii, so its sets are
complete. The probe checks **set inclusion and equality**, not just these
cardinalities, against depth 4. Therefore the earlier equal sets are complete
as well; smaller sets at every preceding depth establish minimal catalog
depth. This reasoning upgrades the observed counts to certified values.

| Collar radius n | Minimum coverage depth, proved | Minimum catalog depth, proved | Legal collar count |
|---:|---:|---:|---:|
| 1 | 1 | 1 | 14 |
| 2 | 2 | 2 | 29 |
| 3 | 3 | 2 | 50 |
| 4 | 3 | 3 | 77 |
| 5 | 4 | 3 | 110 |
| 6 | 4 | 3 | 149 |

Depth 4 also gives the exact minimum coverage depth for n=7 through 11, since
`R(3)=4` and `R(4)=11`; their catalogs were not counted here. No general formula
for catalog size or catalog depth is inferred from the first six radii.

Growth is not uniform across source vertex types. With supported reconstruction
the minimum/maximum output sizes over the 19 inflated stars are 21/33 at
depth 1, 120/210 at depth 2, 748/1365 at depth 3, and 4983/9183 at depth 4.
The coverage-radius ranges in the previous table also vary with context.
The guarantee uses their minimum after maximizing over eligible stars at each
target, rather than assuming a common star shape or a uniformly centered collar.

## 8. Reproduction, verification, and limits

`collar_probe.cpp` is a standalone research utility, not a production collar
enumerator. Its small `partial` evaluator models the information content
assumed in the question; it does not change or repair the production converter.
No application command, API, renderer, or CMake target was changed.

Build from an x64 MSVC developer shell at the repository root:

```text
cl /nologo /std:c++latest /EHsc /O2 /MD /Isrc research\collar_probe.cpp src\triangle_grid.cpp src\combinatorial_hex.cpp src\hex_grid.cpp src\finite_patch.cpp /Fe:out\build\collar_probe.exe /Fo:out\build\
out\build\collar_probe.exe extendibility 4
out\build\collar_probe.exe selfcheck 4
out\build\collar_probe.exe admissible 4
out\build\collar_probe.exe contexts 4
out\build\collar_probe.exe geometry 4 4
out\build\collar_probe.exe catalog 4 4
```

For `geometry` and `catalog`, the first number is the depth used to obtain
the certified finite seed data, and the second is the maximum inflation depth
to check. Geometry enumerates output centers in each central parent polygon,
not merely centers well inside an arbitrarily large sampled supertile.

Verified research checks:

* All 132 single-sector cases (22 labels times six sectors) uniquely decode.
* All 57 combinations of 19 vertex types and depths 1–3 agree between direct
  level-dependent placement and repeated repository finite-patch substitution,
  leaf by leaf after a common translation.
* For those cases, polygon footprints and actual substituted triangle
  footprints agree at every output coordinate.
* The finite upper/lower Stage-1 and radius-2 certificates run successfully.
* Coverage is exhausted through depth 4, including explicit failures one
  radius beyond each minimum.
* Catalogs through depth 4 have the set equalities used above.

An independent read-only review reran `selfcheck`, `admissible`, `geometry`, and
`catalog` with the arguments above, reproduced the results, and found no
correctness issues in the finite reductions or the extension argument.

Saved outputs are `collar-certificate.txt`, `collar-contexts.txt`,
`collar-extendibility.txt`, `collar-selfcheck.txt`, `collar-geometry.txt`, and
`collar-catalogs.txt`, all beside this report.

The available CTest run passed four of five suites. `hex_grid_reconstruction`
reached its final SVG cleanup and failed removing the temporary SVG on Windows;
the test keeps its `std::ifstream` open at the `filesystem::remove` call.
The other four suites were rerun and passed. This unrelated test cleanup was
left unchanged. The research program was separately built against the current
exact sources with optimization enabled.

The remaining limitations are explicit: no closed-form exact R(k) for arbitrary
k, no production implementation of the partial converter or collar enumerator,
and no theorem that local E/V compatibility alone implies global legality.
Neither additional edge seeds nor larger enumeration seeds are needed to
restore completeness. The necessary addition to a naive implementation is
the **certified Stage-1 catalog and the coverage calculation**, with correct
partial reconstruction and phase registration. Larger source radius-2 patches
are used only for the proof/verification of coverage, not as Stage-2 seeds.
