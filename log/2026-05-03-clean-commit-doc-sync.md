# Log: Clean Commit Documentation Synchronization

Date: 2026-05-03

## Summary

Completed a documentation cleanup pass after auditing the repository state for a
clean commit.

## Intended Scope

- Synchronize stale internal vignette text with the current package API.
- Re-render internal vignette HTML outputs from the updated sources.
- Keep untracked local scratch artifacts out of the commit scope through local
  git excludes.
- Validate with R source parsing, Rd checks, and testthat.

## Changes

- Updated `internal/vignettes/02-matern-workflow.Rmd` and its `.rmarkdown`
  sibling to describe log-link `backend = "auto"` as Fisher Laplace, not the
  old sparse observed-Hessian Laplace default.
- Updated the Matern workflow PC-prior section to reflect current behavior:
  identity-link `backend = "auto"` remains on the exact Gaussian backend while
  explicit `backend = "inla"` is the independent INLA/SPDE comparison route.
- Updated fixed-hyperparameter documentation to state that `fix_g = TRUE`
  fixes smoothing hyperparameters from `g_init`, while `beta_fixed` fixes the
  intercept separately.
- Updated `internal/vignettes/03-package-sanity-check.Rmd` so PC-prior sanity
  checks assert stable current fields rather than removed Step A/Step B
  diagnostic equalities.
- Re-rendered all internal vignette HTML outputs in
  `internal/vignettes/rendered/`.
- Added `.vscode/`, `attempts_in_EBMF/`, and `tryingSpatial.R` to local
  `.git/info/exclude` so scratch artifacts stay out of normal commit status
  without deleting local files.

## Notes

- Core package source, public API, tests, and package Rd files were not targeted
  for behavioral changes.
- The cleanup intentionally does not commit or remove local experiment files.
- A full `testthat::test_dir("EBSmoothr/tests/testthat")` run was started, but
  was interrupted after the INLA external binary entered the known repeated
  Hessian diagnostic loop. Targeted test files were run separately afterward
  with output redirected to `/tmp` logs.

## Validation

- `Rscript internal/vignettes/render_internal_vignettes.R`
  - Passed; regenerated all four rendered internal vignette HTML files.
- R source parse and Rd check:
  - Passed for all `EBSmoothr/R/*.R` and `EBSmoothr/man/*.Rd` files.
- Targeted testthat files:
  - `test-eb-smoother.R`: passed.
  - `test-lgp.R`: passed.
  - `test-spatial-scores.R`: passed.
  - `test-matern.R`: passed; INLA emitted its known external binary retry /
    abort-trap diagnostics, but testthat completed successfully.
