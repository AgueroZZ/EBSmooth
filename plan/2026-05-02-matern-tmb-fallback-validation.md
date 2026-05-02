# Matern TMB Repair Validation Plan

Date: 2026-05-02

## Goal

Prevent log-link Matern fits from returning non-finite posterior moments while
keeping the default Laplace backend on the TMB path whenever possible.

## Steps

1. Reproduce the bad mixed-backfit update and save the finite `x` and `s`
   pseudo-data that trigger non-finite posterior moments.
2. Compare the default/TMB Laplace backend with the R reference Laplace backend
   on the same pseudo-data to separate input-data problems from TMB optimizer
   problems.
3. Treat TMB fits with non-finite log-likelihoods, non-finite Step B objectives,
   or non-finite posterior moments as numerically invalid.
4. Improve the TMB path before considering R fallback:
   - add a larger-range Step A start for known-noise empirical-Bayes log-link
     fits;
   - refine the random-effect mode with Step B when Step A reports nonzero
     convergence but has a usable mode;
   - test, but do not keep, a TMB-only fixed-`g` repair that slightly shrinks
     the fitted spatial marginal scale when reference comparison rejects it.
5. Keep R fallback only as a last-resort default-backend safety net. Explicit
   `backend = "laplace_tmb"` should either return a valid TMB fit or report the
   validation failure.
6. Add focused validation tests for the TMB-fit validation helper and verify the
   captured bad pseudo-data agrees between `laplace_tmb` and `laplace_r` when no
   automatic repair is used.
7. Compare any proposed TMB repair against `laplace_r`; reject finite-only
   repairs when they produce posterior moments that are not reference-consistent.
