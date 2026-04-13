## Update

Added a lightweight internal vignette to sanity-check the current public
package interface.

## Completed Changes

- Added `internal/vignettes/03-package-sanity-check.Rmd`.
- Covered the following public fitting paths:
  - L-GP with identity link;
  - L-GP with log link;
  - L-GP with fixed smoothing parameter and fixed beta coefficients;
  - exact Matern in one dimension;
  - exact Matern in two dimensions;
  - Matern with `pc.penalty` in Step A mode;
  - Matern with `pc.penalty` and fixed hyperparameters.
- Added simple assertions so the vignette fails to render if key invariants are
  broken.
- Updated `internal/vignettes/README.md`.

## Validation

- `Rscript internal/vignettes/render_internal_vignettes.R`

## Notes

- This vignette is intentionally light and complements, rather than replaces,
  the formal `testthat` suite.
