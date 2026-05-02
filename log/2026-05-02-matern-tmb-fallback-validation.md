# Matern TMB Repair Validation Log

Date: 2026-05-02

## Diagnosis

The simulation mixed backfit failed after a log-link Matern loading update
returned `Inf` posterior moments. The pseudo-data passed into the Matern EBNM
were finite (`x` ranged roughly from `-0.42` to `1.42`, with constant
`s = 0.1676`), so the failure was not caused by missing or infinite input data.

On the captured pseudo-data, the original TMB Laplace path could return `false
convergence (8)`, a non-finite Step B joint objective, and non-finite posterior
moments. The inputs were not the problem: the R reference Laplace path returned
a finite fit. Additional TMB-only checks showed that a better Step A start can
find the high-likelihood region, but the unpenalized log-link objective can sit
near a singular conditional Hessian where the lognormal posterior moments
overflow.

## Changes

- Added `.matern_laplace_fit_has_finite_posterior()`.
- Added `.matern_laplace_tmb_invalid_reason()`.
- Added an extra larger-range Step A start for known-noise empirical-Bayes
  log-link TMB fits, even when `matern_n_starts = 1`.
- Updated Step B handling so a non-finite Step B objective does not immediately
  discard a usable Step A mode. The fit is marked invalid and can proceed to
  TMB repair.
- Evaluated a TMB-only fixed-`g` repair for known-noise fits. The repair slightly
  shrank the fitted spatial marginal scale and returned the first TMB candidate
  with finite posterior moments.
- Follow-up validation showed that the fixed-`g` sigma-shrink repair was not
  reliable: it could turn `Inf` posterior moments into finite but astronomically
  large posterior moments that did not agree with the R reference implementation.
  The repair code was removed from the automatic path.
- With the repair disabled, the captured bad pseudo-data gives matching
  `laplace_tmb` and `laplace_r` results for log likelihood, Matern
  hyperparameters, beta, and posterior moments. This points to instability in
  the unpenalized log-link Laplace posterior moments rather than a simple TMB
  implementation mismatch.
- Default `backend = "laplace"` still has an R reference fallback as a last
  resort if TMB cannot produce a valid fit.
- Explicit `backend = "laplace_tmb"` now reports validation failures for
  non-finite TMB results instead of silently returning them.
- Added a focused test for the TMB validation helper.

## Validation

- With automatic repair disabled, the captured bad pseudo-data returns matching
  `laplace_tmb` and `laplace_r` values: log likelihood about `107.78`, Matern
  range about `0.471`, Matern sigma about `1.707`, and beta about `-4.824`.
- The posterior moments are finite but extremely large under both
  implementations, so finite posterior moments alone are not a sufficient
  reliability criterion for this unpenalized log-link model.
- Parse checks passed for `EBSmoothr/R/02_Matern.R` and
  `EBSmoothr/tests/testthat/test-matern.R`.
- Direct validation-helper assertions passed.
