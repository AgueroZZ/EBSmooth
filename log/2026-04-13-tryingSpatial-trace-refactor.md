## Update

Refactored `tryingSpatial.R` into reusable analysis helpers and added a local
trace recorder for smooth-EBNM hyperparameter dynamics inside `flashier`.

## Completed Changes

- Replaced the one-off script structure with helper functions for:
  - simulating the spatial toy dataset;
  - building spatial and non-spatial `ebnm_list` objects;
  - running `flashier` with consistent convergence and verbose settings;
  - normalizing and aligning estimated components to truth;
  - plotting line recoveries and interpolated factor surfaces.
- Added `make_flash_trace_recorder()` to trace `flashier:::solve.ebnm` and
  collect per-update records containing:
  - `phase`, `factor_k`, and side;
  - `rho`, `sigma`, `beta`, and their warm-start deltas;
  - block log-likelihood and KL;
  - posterior-mean change magnitude.
- Added `summarize_matern_trace()` and `plot_matern_trace_grid()` to make the
  trace directly inspectable from the analysis script.
- Kept the spatial and standard-flash comparisons as top-level analysis steps,
  but routed them through the new helpers.

## Validation

- `Rscript -e 'parse(file = "tryingSpatial.R")'`
- Ran a reduced toy example by evaluating only the helper definitions and then
  calling `run_flash_model(..., trace_matern = TRUE)` with:
  - `grid_side = 4`;
  - `n = 20`;
  - `greedy_Kmax = 2`;
  - `backfit_maxiter = 2`.
- Confirmed that the trace summary contains nontrivial `rho`, `sigma`, `beta`,
  and posterior-change trajectories for the Matern updates.

## Notes

- The local validation still prints pre-existing `sh: /bin/kstat: No such file
  or directory` messages from the runtime environment; these did not block the
  trace collection.
- During greedy fitting, the temporary factor index may not be fully populated
  until the factor is inserted into the `flash` object, so early trace rows
  should be interpreted together with `phase` and `trace_step`.
