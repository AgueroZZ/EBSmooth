# Matern PC Prior Sensitivity Check

## What changed

- Added a fixed-data sensitivity script that compares the exact EB fit to several legacy INLA fits with different joint PC-prior strengths.
- Used a slightly smoother and noisier 2D dataset so the check remains visible but lightweight.

## Why

The earlier exact-vs-legacy comparison showed near-agreement under one reasonable PC-prior setup. This sensitivity check asks a narrower question: how much separation appears once the PC prior is intentionally weakened or strengthened on the same data?

## Key findings

- On the fixed dataset, the exact EB baseline fit was `(range, sigma, beta0) = (0.2914, 0.3527, 0.1959)`.
- Under the weakest balanced PC prior, the legacy INLA fit moved to `(0.3035, 0.3630, 0.1958)` with exact-objective gap `0.0277`.
- Under the strong smooth/low-variance prior, the legacy INLA fit moved to `(0.3107, 0.3641, 0.1956)` with exact-objective gap `0.0585`.
- Under the strong rough/high-variance prior, the legacy INLA fit moved to `(0.3006, 0.3626, 0.1959)` with exact-objective gap `0.0242`.
- Even in the strongest tested scenarios, the legacy posterior means remained almost indistinguishable from the exact EB posterior mean: correlations above `0.99996` and RMSE below `0.0028`.
- For these tested priors, the legacy INLA fit did not materially separate from the exact EB fit on this dataset.
