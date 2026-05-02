## Objective

Refactor `tryingSpatial.R` into reusable helpers and add local tracing for
smooth-EBNM hyperparameter dynamics during `flashier` updates.

## Plan

1. Replace the monolithic analysis script with helper functions for:
   - simulation;
   - flash fitting;
   - component normalization and alignment;
   - line and surface plotting.
2. Add a local trace recorder that hooks into `flashier:::solve.ebnm` and
   stores per-update diagnostics such as:
   - factor index and side;
   - phase (`greedy` or `backfit`);
   - `rho`, `sigma`, and `beta`;
   - block log-likelihood and KL;
   - posterior-mean change magnitude.
3. Expose the trace as a regular data frame and add summary/plot helpers for
   inspecting plateau behavior.
4. Validate on a smaller toy example before leaving the full interactive script
   in place.
