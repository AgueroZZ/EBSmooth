# Matern Flat-Beta and Warm-Start Fix Log

## Changes
- Corrected Matern TMB flat-prior beta handling: `beta_prec = 0` now includes beta in the TMB random-effect block and adds no beta prior term.
- Kept empirical-Bayes beta as an optimized fixed parameter and proper-prior beta as an integrated random effect with a Gaussian prior.
- Reused `objA$env$last.par.best` from TMB Step A as the final random-effect mode, avoiding the slow cold-start Step B in normal cases.
- Added Step B fallback diagnostics through `tmb_mode_source`.
- Warm-started log-link INLA comparable Laplace objective evaluations from the INLA spatial-field mode.
- Updated tests to require TMB support for flat-prior beta and to check the INLA warm-start likelihood path.

## Validation Notes
- Flat-prior beta was manually checked before implementation by comparing a direct TMB random-beta objective against the R Laplace flat-beta reference.
- Runtime validation should be repeated after focused tests because the previous `n ~= 1000` TMB timing was dominated by the removed cold-start Step B.
- Focused validation passed: `test-matern.R`, `test-eb-smoother.R`, `test-lgp.R`, `R CMD INSTALL EBSmoothr`, and Rd checks for the updated Matern docs.
- A local runtime smoke after the fix gave 1D `n = 1000` TMB runtimes of about 0.47 seconds for fixed beta and 0.30 seconds for EB beta, with log likelihoods aligned to INLA after warm-start correction.
- The same smoke showed 2D `n = 1024` TMB runtimes around 22-28 seconds. A timing breakdown for 2D fixed beta found the main remaining cost in TMB Step A sparse AD/logdet optimization, not posterior summaries.
