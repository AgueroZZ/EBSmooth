# Simplify EBSmoothr Matern Backend API

## Summary
Simplify the public Matern backend interface so users choose high-level
backends while older low-level names remain compatibility aliases.

## Implementation
- Document public Matern backends as `auto`, `exact`, `laplace`, `laplace_r`,
  and `inla`.
- Map `laplace_tmb` to `laplace`; the actual implementation is reported in
  `laplace_implementation`.
- Map `inla_pc` to `inla` when `pc.penalty` is supplied; keep the existing
  error when `inla_pc` is requested without `pc.penalty`.
- Report `fit$backend` as `exact`, `laplace`, or `inla`.

## Validation
- Parse EBSmoothr source and targeted tests.
- Run targeted `test-eb-smoother.R` and `test-matern.R`.
- Parse SpatialEBMF `code/simulation_functions.R`.
