# Softplus Reference Rename and Runtime Validation

## Summary

Completed a follow-up cleanup after the softplus posterior moment fix. Internal
constant-baseline helpers no longer use the old nonspatial naming, and spatial
selection outputs now use reference-route labels for non-Matern alternatives.

## Changes

- Renamed internal `eb_smoother()` constant-baseline helpers from
  `.nonspatial_*` to `.constant_*`.
- Updated `spatial_select()` score and permutation outputs so the non-spatial
  route is reported as `selection = "reference"` and
  `route_label = "reference"`.
- Added `loglik_reference` and `reference_fits` outputs.
- Kept `loglik_nonspatial`, `nonspatial_beta`, `nonspatial_noise_sd`, and
  `nonspatial_fits` as compatibility aliases for older scripts.
- Regenerated Rd documentation for `spatial_scores()` and `spatial_select()`.

## 2D Matern Runtime Audit

Runtime was checked through the public `eb_smoother()` path on fixed-noise 2D
Matern fits with TMB Laplace backends.

For an 8 by 8 grid (`n = 64`) with a reused 2D Matern setup (`132` mesh
vertices):

- `link = "softplus"`, `backend = "auto"` resolved to observed Laplace and had
  median elapsed time `0.784` seconds across three repeats.
- `link = "softplus"`, `backend = "laplace_fisher"` had median elapsed time
  `0.762` seconds across three repeats.
- `link = "log"`, `backend = "auto"` resolved to Fisher Laplace and had median
  elapsed time `1.599` seconds across three repeats.
- `link = "log"`, `backend = "laplace"` had median elapsed time `1.588`
  seconds across three repeats.

For direct end-to-end calls using `locations` and mesh construction inside
`eb_smoother()` on the same grid:

- Softplus auto elapsed time was `0.991` seconds.
- Log auto elapsed time was `2.098` seconds.

For a small scaling smoke with reused setup:

- `n = 64`: softplus auto `1.036` seconds; log auto `1.587` seconds.
- `n = 100`: softplus auto `1.426` seconds; log auto `2.916` seconds.

The observed timings are small-problem smoke checks, not a comprehensive
benchmark suite. They indicate that the softplus 2D Matern path is working and
has reasonable runtime relative to the log-link Matern path.

## Validation

- `roxygen2::roxygenize("EBSmoothr")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-softplus-moments.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-lgp.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-eb-smoother.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-spatial-scores.R")`
- `R CMD build EBSmoothr --no-build-vignettes`
- `R CMD check EBSmoothr_0.1.0.tar.gz --no-manual --no-build-vignettes --no-tests`

The targeted test files completed with zero failures. The source package check
completed with the existing `unlockBinding()` NOTE in `R/02_Matern.R` and no
errors or warnings.
