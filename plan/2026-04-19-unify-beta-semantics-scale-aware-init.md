# Unified Beta Semantics And Scale-Aware Initialization

## Goal
- Unify the public `beta` semantics across `eb_smoother()`, `ebnm_Matern_generator()`, and `ebnm_LGP_generator()`.
- Make default Matérn initialization data/prior-scale-aware when `g_init = NULL`.
- Restrict empirical-Bayes `beta` for Gaussian Matérn + PC prior to the exact backend.

## Planned Changes
- Add public `beta_prec` support to `eb_smoother()`.
- Extend `Matern`, `LGP`, and `nonspatial` family objects to carry beta state.
- Add `Nonspatial(beta, beta_prec)` family object.
- Make Matérn exact paths support:
  - fixed beta
  - empirical-Bayes beta
  - flat beta prior
  - proper beta prior
- Make Matérn `backend = "auto"` choose:
  - `exact` when `pc.penalty = NULL`
  - `exact` when `pc.penalty != NULL` and beta is empirical-Bayes
  - `inla_pc` when `pc.penalty != NULL` and beta is fixed or prior-based
- Reject `backend = "inla_pc"` with empirical-Bayes beta.
- Translate `LGP_setup$betaprec` into the new public semantics as a legacy fallback.
- Update tests to target the new semantics instead of the previous `profile` / `integrated_flat` naming.

## Verification Plan
- Load the package with `pkgload::load_all("EBSmoothr")`.
- Run focused tests for:
  - Matérn beta modes and backend dispatch
  - `eb_smoother()` wrapper behavior
  - LGP beta semantics and legacy fallback
