## Objective

Run internal simulation validation for the exact Matern backend in
`EBSmoothr`, focusing on three questions:

1. whether the empirical-Bayes smoothing hyperparameters recover the truth in a
   well-specified simulation setting;
2. whether the fitted posterior surface tracks the true simulated surface in
   one and two dimensions;
3. whether the exact marginal likelihood implementation is numerically correct.

## Plan

1. Build an internal validation script that simulates data directly from the
   package's exact Matern latent-field model.
2. Run a one-dimensional parameter-recovery experiment over increasing sample
   sizes and summarize the fitted `range`, `sigma`, and surface error metrics.
3. Run one-dimensional and two-dimensional visual surface-recovery examples and
   save figures.
4. Check the exact marginal likelihood against a dense Gaussian evaluation on a
   small problem.
5. Save results under `internal/simulations/` and record the outcome in the
   collaborator log.
