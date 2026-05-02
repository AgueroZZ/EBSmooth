# Runtime Benchmark: 2D Exact Matern vs INLA-PC

- Date: 2026-04-17
- Data geometry: 60 x 50 regular grid (3000 observations)
- Matern setup max.edge: 0.1, 0.2
- SPDE nodes: 3352

## Results

- Exact known-noise fit: 7.482 s (completed)
- Exact learned-noise fit: 9.198 s (completed)
- INLA-PC learned-noise fit: 9.147 s (completed)
- Single exact objective evaluation: 0.040 s (completed)

## Notes

- The exact objective evaluation isolates the shared sparse-factorization path used by the 2D exact optimizer.
- Compare this file against earlier slow-path timings to catch regressions in permutation-aware sparse linear algebra.
