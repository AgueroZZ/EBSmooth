# Exact Gaussian Learned-Noise Scale Extension

## Summary

Implemented internal support for exact-Gaussian learned-noise Matern fits with
fixed per-observation noise scales. This prepares the Gaussian engine needed
for a future Fisher-PQL backend where pseudo-response variance is
`noise_sd^2 / g_i^2`.

## Changes

- Added an internal effective-noise helper with a preserved
  `noise_scale = NULL` fast path.
- Threaded `noise_scale` through exact learned-noise objective helpers and the
  fitted-object construction.
- Stored both scalar `fitted_noise_sd` and effective `fitted_s` on exact
  learned-noise fits.
- Added focused tests for unit-scale equivalence, heterogeneous noise scales,
  objective equivalence, and invalid scale validation.
- Updated `NEWS.md` with the internal infrastructure and validation note.

## Runtime Notes

- Pre-change local smoke timing for the default exact learned-noise path
  (`noise_scale = NULL`, `n = 80`) was `0.276, 0.238, 0.227` seconds across
  three warmed runs.
- Post-change local smoke timing for the same default path was
  `0.274, 0.231, 0.239` seconds across three warmed runs. The median changed
  from `0.238` seconds to `0.239` seconds, within the 5% overhead guard.

## Validation

- Parsed `EBSmoothr/R/02_Matern.R`, `EBSmoothr/R/03_eb_smoother.R`, and
  `EBSmoothr/tests/testthat/test-matern.R`.
- Ran a direct objective-equivalence smoke check for heterogeneous
  `noise_scale` across empirical-Bayes, fixed, flat-prior, and proper-prior
  beta modes.
- Ran `testthat::test_dir("EBSmoothr/tests/testthat", filter = "matern")`:
  `298` passed, `0` failed, `0` skipped. Local INLA emitted optimizer fallback
  and subprocess crash messages during existing inlabru/INLA cases, but
  testthat completed successfully.
