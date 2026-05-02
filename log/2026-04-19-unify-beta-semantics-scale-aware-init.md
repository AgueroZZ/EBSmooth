# Update Log: Unified Beta Semantics And Scale-Aware Initialization

## Code Changes
- Added `beta_prec` handling to `eb_smoother()`.
- Added `Nonspatial(beta, beta_prec)` and exported it through `NAMESPACE`.
- Extended `Matern(theta, sigma, beta, beta_prec)` and `LGP(scale, beta, beta_prec)` state handling.
- Implemented exact Matérn support for:
  - empirical-Bayes beta
  - fixed beta
  - flat beta prior
  - proper beta prior
- Added scale-aware Matérn initialization helpers for:
  - `theta_init`
  - `sigma_init`
  - `beta_init`
  - `noise_sd_init`
- Updated Matérn wrapper/generator dispatch so empirical-Bayes beta with PC prior uses the exact backend, while fixed/prior beta can use `inla_pc`.
- Updated `matern_objective_breakdown()` to recognize the new public beta mode names.
- Updated nonspatial exact fits so fixed / EB / flat / proper beta modes all return `fitted_g` with beta state.
- Updated LGP generator and learned-noise wrapper to use:
  - `beta_prec = NULL` for empirical-Bayes beta
  - `beta_prec = 0` for flat beta prior
  - `beta_prec > 0` for proper beta prior
  - `beta_fixed` even when `fix_g = FALSE`

## Test Changes
- Rewrote `test-matern.R`, `test-eb-smoother.R`, and `test-lgp.R` around the new beta semantics and backend rules.

## Notes
- Focused package loading succeeds with `pkgload::load_all("EBSmoothr")`.
- The test suite was updated to the new semantics; long-running INLA-based checks remain noisy on this machine because of repeated `/bin/kstat` warnings from upstream dependencies.
