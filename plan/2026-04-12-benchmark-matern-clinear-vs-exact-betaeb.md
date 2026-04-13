## Objective

Benchmark the computational efficiency of the `clinear`-based INLA prototype
for the Matern `betaprec < 0` case against the current exact profiled
implementation.

## Plan

1. Build a reproducible benchmark with a shared spatial setup for both methods.
2. Compare fit-phase runtime rather than setup time so the benchmark matches
   the package's generator/fitter usage pattern.
3. Use one moderate 1D example and one moderate 2D example.
4. Record runtime together with objective and parameter gaps to confirm the two
   methods still solve nearly the same optimization problem.
5. Save the benchmark outputs under `internal/simulations/results/`.
