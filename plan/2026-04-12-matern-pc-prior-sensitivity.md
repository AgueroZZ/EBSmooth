# Matern PC Prior Sensitivity Check

## Goal

Check, on a fixed dataset, whether the legacy INLA fit begins to separate from the exact EB fit when the PC prior is made weaker or stronger.

## Planned Changes

- Simulate one fixed 2D dataset from the exact Matern model.
- Fit the current exact EB Matern model once.
- Refit the legacy INLA model several times while varying the PC-prior tail probabilities.
- Compare fitted hyperparameters, runtime, posterior-surface similarity, and the exact objective evaluated at the legacy fitted hyperparameters.

## Deliverables

- New script `internal/simulations/study_matern_pc_prior_sensitivity.R`
- CSV summaries and a small figure under `internal/simulations/results/`
- A markdown/HTML note summarizing whether the fits separate materially under stronger priors
