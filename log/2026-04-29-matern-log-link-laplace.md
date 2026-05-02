# Matérn Log-Link Sparse Laplace Extension

Implemented known-noise positive Matérn smoothing with `link = "log"`.

## Changes
- Added a package-owned sparse Laplace backend for Matérn known-noise fits.
- Added an INLA backend that supports unpenalized and PC-prior Matérn fits.
- Updated `ebnm_Matern_generator()` and `eb_smoother()` backend dispatch:
  `backend = "auto"` keeps exact Gaussian identity-link behavior and uses
  manual Laplace for log-link Matérn fits.
- Kept learned-noise Matérn log-link fits unsupported with an explicit error.
- Added tests for manual/INLA log-link agreement, beta modes, wrapper behavior,
  and identity-link Laplace regression against the exact backend.

## Notes
- For INLA log-link fits, the stored primary `log_likelihood` is recomputed
  with the package Laplace objective at the INLA mode so it is directly
  comparable to the manual backend. Raw INLA objective fields remain stored as
  diagnostics.

## Follow-up: Backend Equivalence and Runtime
- Found that `backend = "inla"` was not equivalent to the default log-link
  manual Laplace backend under default beta handling. The manual path uses
  empirical-Bayes beta optimization when `beta_prec = NULL`, while INLA treats
  an unfixed fixed effect as part of the latent Gaussian model with a flat
  prior. These are different objectives under the log link.
- Initially updated Matern INLA dispatch to reject empirical-Bayes beta mode
  explicitly, so different beta conventions could not be silently compared.
  This temporary guardrail was superseded by the unified beta semantics
  follow-up below.
- Added a flat-prior beta agreement test for log-link manual Laplace versus
  INLA, in addition to the existing fixed-beta comparison.
- Reduced manual Laplace overhead by skipping posterior variance and
  latent-field summary computations during outer hyperparameter optimization.
  Those quantities are now computed only once at the final fitted mode.
- Tried a sparse Newton inner solver for the log-link latent mode, but did not
  keep it because the non-convex likelihood can lead Newton to a different
  latent mode than the robust `nlminb` path.

## Validation
- `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")` passed.
- In a 60-point log-link empirical-Bayes beta example, manual Laplace runtime
  was about 7.8 seconds after skipping repeated posterior summaries; explicit
  `backend = "inla"` now errors with guidance to use `beta_fixed` or
  `beta_prec = 0`.
- In a 120-point log-link flat-prior beta example, manual Laplace and INLA
  agreed on the comparable Laplace objective and posterior means:
  manual `log_likelihood = 37.58064`, INLA recomputed
  `log_likelihood = 37.58062`, max posterior-mean difference about `0.00255`.
  Runtime was about 17.3 seconds for manual Laplace versus about 4.1 seconds for
  INLA on this machine.

## Follow-up: Unified Beta Semantics Across Backends
- Removed the temporary empirical-Bayes beta rejection for `backend = "inla"`
  and `backend = "inla_pc"`.
- Added INLA external beta profiling helpers. Empirical-Bayes beta now runs a
  small local profile over fixed-beta offset INLA fits, then returns the final
  optimized fixed-beta INLA fit as an empirical-Bayes beta result.
- Applied the same profiling approach to known-noise INLA fits and learned-noise
  `inla_pc` fits.
- Preserved fixed-beta and beta-prior INLA behavior unchanged.
- Updated documentation to describe INLA empirical-Bayes beta as external beta
  profiling.
- Updated tests so INLA empirical-Bayes beta is expected to work, and added a
  log-link empirical-Bayes beta manual-vs-INLA agreement check.
- Validation:
  `testthat::test_file("EBSmoothr/tests/testthat/test-matern.R")` passed, and
  `tools::checkRd()` passed for the updated Matern and `eb_smoother` docs.
