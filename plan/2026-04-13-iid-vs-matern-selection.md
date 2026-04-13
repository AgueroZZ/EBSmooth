## Objective

Implement a reproducible internal study that compares exact Matern smoothing
against an iid normal prior under iid-normal truth, and quantify when raw
optimized marginal likelihood selects the iid model.

## Plan

1. Add a dedicated simulation script for the iid-vs-Matern selection study.
2. Fit two candidate models on every replicate:
   - exact Matern with optimized `(beta0, range, sigma)`;
   - `ebnm_normal(mode = "estimate", scale = "estimate")`.
3. Use the raw optimized marginal likelihood difference
   `delta = loglik_Matern - loglik_iid` with tolerance `1e-6` to label each
   replicate as `iid`, `matern`, or `tie`.
4. Cover both fixed-domain and increasing-domain regimes across the requested
   `n` and `tau0` grids.
5. Run the study in two stages:
   - a pilot pass over every design cell;
   - a confirmatory pass for the baseline cell and pilot cells whose iid
     selection frequency lies in `[0.35, 0.65]`.
6. Save raw draws, cell-level summaries, a markdown summary, and the requested
   figures under `internal/simulations/results/`.
7. Add a mesh-sensitivity follow-up for Matern-selected replicates with
   fitted range below the mesh edge.
8. Validate the implementation with a smoke run before leaving the branch.
