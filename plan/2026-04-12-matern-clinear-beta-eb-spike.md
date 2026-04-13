## Objective

Test whether an INLA `pcmatern` model with the intercept represented through
`clinear` can recover the same profiled objective and fitted parameters as the
current exact Gaussian Matern implementation when the intercept is treated as
an empirical-Bayes parameter.

## Plan

1. Build a small 1D and 2D prototype using `inla.spde2.pcmatern()` for
   `(range, sigma)` and `clinear` for the intercept.
2. Verify the correct `clinear` encoding for a pure intercept term.
3. Compare the INLA Step A penalized objective against the exact profiled
   Gaussian objective plus the PC-prior contribution.
4. Compare the fitted `(beta0, range, sigma)` values against the exact optimum.
5. Save the spike as an internal simulation record for future reference.
