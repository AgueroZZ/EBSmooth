# Non-Identity Initialization Scale Fix

## Summary

Implemented link-consistent initialization for Matern and L-GP log/softplus
fits. Default beta starts now match response-scale means; latent prior scale
starts are kept on eta scale; learned common observation noise starts are kept
on the raw response scale.

## Versioning

- Recorded this fix under `EBSmoothr 0.2.4`. `EBSmoothr/DESCRIPTION` and
  `EBSmoothr/NEWS.md` already carry the `0.2.4` package version for the current
  minor update, so no additional version bump was made.

## Changes

- Added shared internal response-scale initialization helpers in
  `EBSmoothr/R/00_setup.R`.
- Updated Matern initialization in `EBSmoothr/R/02_Matern.R` to separate
  latent `sigma_data` from raw learned-noise `noise_sd_init`.
- Added one shared L-GP initialization resolver in `EBSmoothr/R/01_LGP.R` and
  routed the L-GP known-noise, R Laplace/Fisher, and learned-noise paths through
  it.
- Added targeted regression tests in the Matern and L-GP test files.

## Validation

- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet = TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-lgp.R")'`
  passed.
- `Rscript -e 'pkgload::load_all("EBSmoothr", quiet = TRUE); testthat::test_file("EBSmoothr/tests/testthat/test-matern.R", desc = "Matern non-identity initialization keeps latent and raw scales separate")'`
  passed.
- The full `test-matern.R` file was also run with `pkgload::load_all()`. The new
  initialization test passed, but the file still has eight failures in the
  pre-existing `"Matern log-link Fisher backend supports beta modes and learned
  noise"` test because the current dirty worktree routes learned-noise log-link
  Matern fits to `backend = "fisher_pql"` while that test still expects
  `backend = "laplace_fisher"`.
- Plain `testthat::test_file(...)` without `pkgload::load_all()` does not load
  the local package in this checkout and reports missing exported functions; use
  the `pkgload::load_all()` form for local test-file runs.
