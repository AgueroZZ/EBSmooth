# Matern Auto Backend Policy

## Goal

Keep the Matern backends semantically equivalent and use `backend = "auto"` only
for conservative runtime routing. Explicit backend choices must continue to mean
exactly what the caller requested.

## Policy

- Use `exact` for identity-link Matern fits.
- Use package-owned `laplace` for log-link Matern fits by default.
- Use `inla` automatically only for unpenalized 2D log-link fits with
  `beta_fixed` and optimized Matern hyperparameters, where local benchmarks
  showed INLA is faster and the comparable package Laplace objective matches
  TMB/R.
- Keep 2D empirical-Bayes, flat-prior, proper-prior, and PC-prior log-link fits
  on package `laplace` by default until those paths have separate speed evidence.
- Preserve explicit `laplace`, `laplace_tmb`, `laplace_r`, `inla`, and `inla_pc`
  routing.

## Validation

- Add generator tests for the 2D fixed-beta auto-INLA path and compare against
  explicit TMB Laplace.
- Add wrapper tests so `eb_smoother()` mirrors the generator policy.
- Re-run Matern, wrapper, and LGP regression tests.
