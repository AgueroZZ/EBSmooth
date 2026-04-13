## Objective

Add an optional `pc.penalty` mode to the Matern backend so the package can
support either exact empirical-Bayes fitting or an INLA PC-prior workflow
without reintroducing the earlier objective mix-up.

## Plan

1. Extend `ebnm_Matern_generator()` with a generator-level `pc.penalty`
   argument and keep the default exact Gaussian EB implementation unchanged.
2. Route `pc.penalty != NULL` and `fix_g = FALSE` through an INLA Step A fit
   with `inla.spde2.pcmatern()`, and use the Step A penalized objective as the
   main `log_likelihood`.
3. Route `pc.penalty != NULL` and `fix_g = TRUE` through an INLA Step B fit
   with fixed hyperparameters and fixed intercept, and define the main
   `log_likelihood` as Step B plus the PC prior log-density on the internal
   parameter scale.
4. Return explicit backend metadata and backend-specific diagnostic quantities
   so exact and PC-prior fits are distinguishable in downstream analyses.
5. Update tests, package documentation, internal math notes, and internal
   vignettes to reflect the new dual-backend Matern semantics.
6. Run targeted validation covering exact mode, PC-prior Step A mode, and
   fixed-parameter PC-prior mode.
