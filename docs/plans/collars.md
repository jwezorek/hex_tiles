# Collar implementation plan

Specification: user request and research/collar-depth.md.
Work in place, independently, without commits or unrelated refactoring.

1. Add regression tests for supported-sector reconstruction, then extend convert with an explicit mode, exact frame phase inference, shared motif decoding, and observed-sector provenance.
2. Promote the certified 19 seeds and 29 contexts into compiled data. Add exact coverage geometry and deterministic rooted semantic collar canonicalization. Test counts, certificates, and catalog/coverage set equality.
3. Add a separate const-catalog SVG writer and runnable example. Verify XML, counts, geometry bounds, deterministic output and errors; render small and radius-6 sheets.
4. Run CTest and review the full changes. Document API, completeness, provenance and resource limits.

Implementation choices: collar keys keep the central hex at zero; only six proper rotations are compared. The exact core uses checked integer arithmetic. Coverage evaluates every three-face source star in each radius-2 context. Certified catalog depths and geometric coverage certificates have distinct metadata. Explicit hard limits bound depth, radius, and stored catalog entries; exceptions never return partial catalogs.

Progress: research and relevant production sources read. Default conversion behavior must be preserved. Existing Windows SVG cleanup failure is already documented in research.

Progress: supported-sector test first failed on default clipping, then passed with the explicit mode. Catalog and SVG API tests first failed on missing APIs, then passed after implementation. Certified counts, all support/invariance tests, exact compiled source data, 18 edge types, R(0..4), catalog-vs-coverage set equalities and radius-7 coverage selection pass. Release CTest: 9/10; sole failure is the baseline Windows open-input-stream SVG cleanup, reproduced before edits. Radius-1 and radius-6 SVGs rasterized with sharp and visually inspected (full sheets plus radius-6 detail). Reviewer dispatched for final review. Hard limits: radius32, depth5, 1M triangle cells per seed, 1M stored semantic labels, 4M footprint scan candidates. Full API and completeness documentation in COLLARS.md.

Additional manual verification: the public runnable example generated radius 12 using the computed depth-5 coverage certificate, returning 509 collars (no count formula assumed). The depth-5 smoke test completed successfully. Final fresh-context code review found no actionable issues.
Final verification: clean Release rebuild succeeded. Full CTest 9/10 in 32.90 seconds; only baseline hex_grid_reconstruction SVG deletion failure. XML validation also passed for the radius-12 sheet (509 panels, 238721 hexes). git diff --check reports no whitespace errors. No commits or branch/worktree changes.
