## Update

Relocated the internal vignette materials out of `EBSmoothr/inst/` and added
saved rendered outputs.

## Files Added

- `internal/vignettes/render_internal_vignettes.R`
- `internal/vignettes/rendered/` HTML outputs
- `plan/2026-04-11-render-internal-vignettes.md`
- `log/2026-04-11-render-internal-vignettes.md`

## Files Updated

- `internal/vignettes/README.md`
- `internal/vignettes/00-overview.Rmd`
- `internal/vignettes/01-lgp-workflow.Rmd`
- `internal/vignettes/02-matern-workflow.Rmd`

## Notes

- The internal documents are now outside the package installation tree.
- The vignette setup can load `EBSmoothr` from source through `pkgload` when the package is not installed.
- Rendered HTML files are intended to be versioned with the source `.Rmd` files for collaborator review.
- `EBSmoothr/R/01_LGP.R` was updated to use `Matrix::Matrix(..., sparse = TRUE)` when building sparse matrices, because the older `as(..., "dMatrix")` conversion failed during vignette execution under the current Matrix version.
