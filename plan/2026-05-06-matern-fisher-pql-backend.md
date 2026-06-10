# Matern Fisher-PQL Backend

## Goal

Add `backend = "fisher_pql"` for non-identity Matern fits and make it the
default `backend = "auto"` path for log and softplus links. The corrected
backend should accelerate slow log/softplus Matern fits by placing Fisher/PQL
inside the Matern Step A objective: the outer optimizer chooses hyperparameters
and beta/noise where applicable, and the inner latent mode is approximated by a
small fixed number of pseudo-Gaussian exact Matern updates.

## Implementation Plan

- Accept `fisher_pql` in Matern backend resolution for `link = "log"` and
  `link = "softplus"` only.
- Build pseudo responses using
  `z = eta + (x - h(eta)) / h'(eta)`.
- Implement Step A with the exact pseudo-Gaussian Matern machinery instead of
  the TMB random-effect Laplace path.
- For known-noise fits, use `s_eff = s / pmax(h'(eta), g_floor)` inside the
  pseudo-Gaussian observation model.
- For learned-noise fits, call the exact learned-noise solver with
  `noise_scale = 1 / pmax(h'(eta), g_floor)`, so scalar `noise_sd` remains the
  optimized outer noise parameter.
- Keep PC priors on scalar range, sigma, and noise parameters.
- Expose `pql_inner_iter` as a public Matern Fisher-PQL tuning argument in
  `ebnm_Matern_generator()` and `eb_smoother()`, with default `3`.
- Route `backend = "auto"` to `fisher_pql` for log and softplus Matern fits;
  keep identity-link auto routing on the exact Gaussian backend.
- Store the pseudo Step A objective separately from the final score.
- Report primary `log_likelihood` as the original-model Fisher/Laplace score
  evaluated at the final PQL mode. This is not a true re-optimized
  original-model Laplace marginal likelihood.

## Validation Plan

- Test backend dispatch and unsupported-link errors.
- Test pseudo-response algebra for log and softplus links.
- Test that the exact one-step PQL mode still matches an R sparse linear solve
  at fixed hyperparameters when `pql_inner_iter = 1` is requested internally.
- Test that the default Fisher-PQL path reports three pseudo-Gaussian updates
  and that public `pql_inner_iter` overrides are respected.
- Test that Fisher-PQL wrappers no longer use the TMB `model_id = 2` Step A
  path and report `laplace_implementation = "exact_fisher_pql"`.
- Test known-noise `ebnm_Matern_generator()` and learned-noise `eb_smoother()`
  paths for log and softplus links.
- Compare posterior means against the existing Laplace backends on small
  deterministic examples.
- Benchmark representative learned-noise non-identity Matern fits against the
  existing Laplace backend and report speed and drift together.
