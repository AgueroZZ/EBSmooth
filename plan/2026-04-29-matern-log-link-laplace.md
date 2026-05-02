# Matérn Log-Link Sparse Laplace Extension

## Summary
- Add known-noise positive Matérn smoothing with `link = "log"`, targeting
  `x_i ~ N(exp(beta + A_i w), s_i^2)`.
- Use a package-owned sparse Laplace backend as the default log-link backend.
- Add an independent INLA backend for identity and log links, with and without
  PC priors, so manual Laplace fits can be cross-checked.

## Implementation
- Keep the existing exact Gaussian identity-link backend unchanged.
- Add `backend = "laplace"` and `backend = "inla"` to Matérn generator and
  `eb_smoother()` dispatch.
- Use `pc.penalty = NULL` for unpenalized EB and add the existing Matérn PC
  prior contribution only when `pc.penalty` is supplied.
- Keep learned-noise Matérn log-link fits out of scope for this first version.

## Validation
- Check identity-link Laplace against the exact Gaussian backend.
- Check log-link manual Laplace against INLA for posterior means, fitted
  hyperparameters, and the comparable Laplace objective at the INLA mode.
- Cover fixed, empirical-Bayes, flat-prior, and proper-prior beta modes under
  the log link.

## Follow-up: Backend Equivalence and Runtime
- Audit why `backend = "inla"` can differ from the default log-link Laplace
  backend under default beta handling.
- Enforce explicit beta semantics for INLA-backed Matern fits so the package
  does not silently compare empirical-Bayes beta optimization against INLA's
  flat-prior latent beta treatment.
- Remove unnecessary posterior variance and latent-field summary work from
  manual Laplace objective evaluations; compute those quantities only for the
  final fitted mode.
- Add focused tests for the corrected INLA beta guardrail and for manual/INLA
  agreement when both use the same flat-prior beta semantics.

## Follow-up: Unified Beta Semantics Across Backends
- Replace the temporary INLA empirical-Bayes beta guardrail with external beta
  profiling.
- For INLA-backed empirical-Bayes beta fits, optimize beta outside INLA and use
  fixed-beta offset INLA fits for each profile evaluation.
- Keep current fixed-beta and beta-prior INLA paths unchanged.
- Preserve comparable log-link `log_likelihood` by evaluating the package
  Laplace objective at the optimized INLA mode.
