# Inspect INLA Matern Objective Quantities

## What changed

- Added a focused internal script to compare Step A and Step B objective-related quantities in the legacy INLA workflow.
- Included a fixed-theta invariance check to see whether Step B keeps any prior-dependent objective once theta is fixed.

## Why

The exact-vs-legacy comparisons showed that Step B `mlik` behaves like the exact conditional marginal likelihood. This follow-up isolates the different INLA quantities so the role of the PC prior is easier to interpret.

## Key findings

- In Step A, the quantities `mlik`, `misc$log.posterior.mode`, `misc$configs$max.log.posterior`, and `joint.hyper[, "Log posterior density"]` all changed when the PC prior was changed.
- In Step A, `misc$log.posterior.mode` and `misc$configs$max.log.posterior` were numerically equal in the tested cases.
- In Step B, `mlik` matched the exact Gaussian marginal likelihood at the fitted hyperparameters up to numerical rounding.
- When `theta` and `beta` were held fixed in Step B, `mlik` and `max.log.posterior` were exactly invariant across different PC-prior choices, and `log.posterior.mode` was identically zero.
- For the tested workflow, there was no evidence that Step B retains a hidden prior-penalized objective once the hyperparameters are fixed.
