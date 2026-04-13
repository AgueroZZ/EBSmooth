# Inspect INLA Matern Objective Quantities

## Goal

Clarify which INLA output quantities in the legacy Matern workflow correspond to:

- Step A hyperparameter-selection quantities,
- Step B conditional marginal likelihood quantities, and
- prior-dependent versus prior-invariant values.

## Planned Changes

- Run the legacy INLA Step A and Step B workflow on one fixed dataset under several PC-prior settings.
- Record `mlik`, `log.posterior.mode`, `max.log.posterior`, and `joint.hyper` values from Step A and Step B.
- Repeat Step B with the same fixed theta/beta under different priors to check which quantities remain invariant.

## Deliverables

- New script `internal/simulations/inspect_matern_inla_objectives.R`
- CSV tables and a markdown/HTML note under `internal/simulations/results/`
