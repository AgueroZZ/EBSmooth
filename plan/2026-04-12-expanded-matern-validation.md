## Objective

Expand the internal exact-Matern validation to address two follow-up goals:

1. use substantially larger two-dimensional simulation settings so surface
   recovery and hyperparameter recovery are easier to judge;
2. visualize the profile marginal likelihood as a function of `(rho, sigma)`,
   showing the true hyperparameters and the optimized fit.

## Plan

1. Extend the internal simulation script with a larger 2D recovery experiment.
2. Increase the 2D example size used in the saved surface-recovery figure.
3. Add a profile marginal-likelihood surface computation over a `(rho, sigma)`
   grid, profiling out `beta0` at each grid point.
4. Save the new tables and figures and update the internal markdown summary.
5. Re-run the simulation script and confirm the updated outputs render
   correctly.
