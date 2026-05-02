# Component-Specific Spatial Screening for Matrix Factorization

## What Changed
- Reworked the matrix-factorization simulation helpers so component screening now returns both a raw selection label and a conservative update route.
- Added aligned screening summaries to compare true component type against the selected spatial/nonspatial route after permutation matching.
- Added mixed smooth-EBMF support by assigning factor-specific loading-side `ebnm.fn` objects inside a single `flash_backfit()` run.
- Kept the nonspatial iterative route tied to the original loading-side EBNM function, while leaving the factor-side EBNM unchanged.
- Updated selective initialization so only spatial-routed components are replaced by the Matern smoother posterior mean.
- Rewrote `attempts_in_EBMF/simulation.rmd` around a `2 spatial + 2 nonspatial` experiment with four benchmark methods:
  - Regular EBMF
  - Selective Init Only
  - Mixed Smooth-EBMF Warm Start
  - All-Spatial Smooth-EBMF Warm Start
- Added aligned recovery summaries and spatial loading heatmaps across the benchmark methods.

## Validation
- Ran a small helper-level smoke workflow:
  - regular EBMF
  - `eb_smoother` screening
  - selective initialization
  - mixed smooth-EBMF warm start
- Checked that the main Rmd settings produce four fitted regular-EBMF components under the chosen seed.
- Rendered `attempts_in_EBMF/simulation.rmd` successfully with `rmarkdown::render(..., output_format = "md_document")`.
