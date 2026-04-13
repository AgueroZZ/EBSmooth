## Update

Ran a focused spike to test whether the Matern intercept can be treated as an
empirical-Bayes parameter inside INLA by using the `clinear` latent model.

## Completed Changes

- Verified that `clinear` should receive the actual covariate values as its
  first input.
- Verified that using `1:n` was the wrong encoding for an intercept term.
- Built a correct prototype with `beta_cov = 1` for every observation.
- Compared the resulting INLA Step A objective with the exact profiled
  Gaussian objective plus the PC-prior contribution.
- Compared the fitted `(beta0, range, sigma)` values against the exact optimum
  in both 1D and 2D examples.
- Saved the spike outputs under `internal/simulations/results/`.

## Validation

- `Rscript internal/simulations/study_matern_clinear_beta_eb.R`

## Notes

- With the correct `clinear` specification, the INLA Step A objective and the
  exact profiled objective agree to around `1e-5` in the tested examples.
- The fitted `(beta0, range, sigma)` values were also nearly identical to the
  exact optimum.
- This suggests that a fast INLA-based implementation for the `betaprec < 0`
  case is likely feasible.
