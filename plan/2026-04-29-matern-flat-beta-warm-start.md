# Correct Matern TMB Flat-Beta Semantics and INLA Warm Starts

## Summary
- Treat `beta_prec = 0` as flat-prior beta in the Matern TMB backend by integrating beta as a random effect with no beta prior term.
- Fix log-link INLA comparable Laplace objective evaluation by warm-starting the R Laplace inner optimizer at the INLA spatial-field mode.
- Remove the main Matern TMB runtime bottleneck by reusing TMB Step A's optimized random-effect mode instead of cold-starting Step B.

## Implementation Plan
- Allow `prior_flat` in the Matern TMB backend for `alpha = 2` known-noise fits.
- Use TMB random effects as follows:
  - fixed beta: `random = "w"`, mapped beta;
  - empirical-Bayes beta: `random = "w"`, optimized fixed beta;
  - flat-prior beta: `random = c("w", "beta")`, no beta prior term;
  - proper-prior beta: `random = c("w", "beta")`, Gaussian beta prior term.
- Extract `objA$env$last.par.best` after TMB Step A and use it as the final random-effect mode for posterior summaries.
- Keep Step B only as a fallback when the saved TMB mode is unavailable or malformed.
- For INLA log-link fits, pass the INLA spatial-field mode as `initial_mode` to the package Laplace objective; append the fitted beta for integrated-beta modes.

## Validation Plan
- Compare `laplace_tmb` and `laplace_r` for flat-prior beta.
- Verify default `backend = "laplace"` uses TMB for flat-prior beta when supported.
- Check log-link INLA fixed-beta likelihood against TMB/R reference on a larger 1D regression case.
- Run focused Matern, eb_smoother, and L-GP tests, plus installation and Rd checks.
