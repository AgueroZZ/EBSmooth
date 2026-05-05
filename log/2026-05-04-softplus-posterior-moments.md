# Softplus Posterior Moment Audit and Fix

## Summary

Implemented deterministic softplus posterior moments under the marginal
Gaussian Laplace approximation. Softplus response moments now use fixed
Gauss-Hermite quadrature instead of the first-order delta method.

## Changes

- Added shared stable softplus, stable sigmoid, cached Gauss-Hermite quadrature,
  and softplus Gaussian moment helpers.
- Routed LGP, Matern, and constant-reference softplus posterior moments through
  the shared deterministic transform helper.
- Fixed the R Matern Laplace observation terms so `link = "softplus"` uses the
  softplus mean, gradient, and observed/Fisher curvature instead of falling
  through to log-link math.
- Kept identity moments and log-link lognormal moments unchanged.
- Confirmed the TMB softplus path compiles with the AD-safe `logspace_add()`
  implementation.
- Regenerated Rd documentation for LGP, Matern, and `eb_smoother`.

## Validation

- `pkgload::load_all("EBSmoothr", quiet = TRUE)`
- `testthat::test_file("EBSmoothr/tests/testthat/test-softplus-moments.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-lgp.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")`
- `testthat::test_file("EBSmoothr/tests/testthat/test-eb-smoother.R")`
- `R CMD build EBSmoothr --no-build-vignettes`
- `R CMD check EBSmoothr_0.1.0.tar.gz --no-manual --no-build-vignettes --no-tests`

The targeted test files completed with zero failures. `R CMD check` completed
with one existing NOTE for `unlockBinding()` calls in `R/02_Matern.R` and no
errors or warnings.
