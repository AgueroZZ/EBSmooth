## Objective

Rewrite the Matern backend in `EBSmoothr` as an exact empirical-Bayes Gaussian
normal-means smoother, and fix the `fix_g` semantics for both the L-GP and
Matern fitters so users can explicitly provide fixed beta parameters.

## Plan

1. Extend the L-GP fitter API with an explicit `beta_fixed` argument and make
   `fix_g = TRUE` consistently fix both the smoothing parameter and the beta
   coefficients.
2. Replace the current INLA-based Matern Step A / Step B approximation with an
   exact Gaussian marginal-likelihood implementation based on SPDE precision
   matrices and outer optimization over `(range, sigma, beta0)`.
3. Restrict the Matern backend to the identity-link Gaussian normal-means case.
4. Update package documentation, internal math notes, and internal vignettes to
   reflect the new Matern semantics and the new fixed-beta interface.
5. Add a `testthat` test suite covering the new exact Matern likelihood,
   posterior summaries, and the fixed-beta behavior for both fitters.
6. Run targeted R validation for the rewritten code and keep an update log for
   collaborators.
