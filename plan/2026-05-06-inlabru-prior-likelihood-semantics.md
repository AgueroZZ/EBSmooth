# Inlabru prior and likelihood semantics cleanup

## Goal

Make the experimental softplus Matern `backend = "inlabru"` explicit about
prior semantics and comparable in scoring output:

- known-noise inlabru fits must not synthesize a PC prior when `pc.penalty = NULL`;
- learned-noise inlabru fits must require an explicit `pc.penalty`;
- primary `log_likelihood` must be the package Laplace objective evaluated at
  the inlabru fitted hyperparameters, not the raw inlabru/INLA marginal
  likelihood.

## Implementation checklist

1. Replace the inlabru-only default PC-prior helper with an explicit policy
   helper that permits no-PC known-noise fits and rejects no-PC learned-noise
   fits.
2. Recompute inlabru primary likelihoods through the existing manual Laplace
   objective helpers and keep raw inlabru marginal likelihood diagnostics in
   separate fields.
3. Update `ebnm_Matern_generator()` and `eb_smoother()` consistently.
4. Update inlabru tests, documentation, and the sanity-check script.
5. Validate with package load, targeted Matern tests, and a small known-noise
   no-PC smoke comparison.
