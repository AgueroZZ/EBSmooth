# Backend Release-Readiness Summary

Generated: 2026-04-30 00:16:35 EDT

## Accuracy

- Successful fits: 14
- Failed or guarded fits: 1
- Maximum successful objective delta: 0.07222
- Maximum successful posterior-mean delta: 0.01971

Failed or guarded accuracy rows are expected only when explicit INLA cannot be validated against the package reference.

## Runtime

- Successful runtime rows: 16
- Fastest successful row: 1d_log_learned_eb_n100_auto
- Slowest successful row: 1d_log_learned_eb_n100_laplace_r

## Auto Policy

- Identity-link Matern: exact.
- Nonspatial: exact.
- Log-link Matern learned-noise: Laplace/TMB.
- Log-link Matern known-noise 2D fixed beta without PC prior: INLA when validation evidence remains consistent.
- R Laplace: reference/debug implementation, not the auto default when TMB is supported.

Accuracy CSV: internal/simulations/results/backend_release_accuracy.csv
Runtime CSV: internal/simulations/results/backend_release_runtime.csv
