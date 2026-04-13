## Update

Benchmarked the `clinear`-based INLA prototype against the exact profiled
Matern implementation for the `betaprec < 0` interpretation.

## Completed Changes

- Added `internal/simulations/benchmark_matern_clinear_vs_exact_betaeb.R`.
- Benchmarked fit-phase runtime after precomputing the shared spatial setup.
- Included a moderate 1D example and a moderate 2D example.
- Saved both the timing table and the summary table under
  `internal/simulations/results/`.
- Included objective and parameter gaps so runtime comparisons can be checked
  against numerical agreement.

## Validation

- `Rscript internal/simulations/benchmark_matern_clinear_vs_exact_betaeb.R`

## Notes

- The benchmark is intentionally focused on fit-time cost because the spatial
  setup is naturally amortized when a generator is reused.
