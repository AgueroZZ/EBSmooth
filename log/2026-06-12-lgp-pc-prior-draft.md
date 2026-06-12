# 2026-06-12 LGP PC Prior Draft

## Context

This draft explores whether a PC prior on the LGP latent process scale can
break the near-flat flashier EBMF scaling ridge observed in
`EBFPCA/scripts/fit_ebmf.R`.

The implemented draft follows the BayesGP/TMB parameterization for a latent
standard deviation `sigma_u = exp(-theta / 2)`:

```r
lambda <- -log(alpha) / anchor
log_prior <- log(lambda / 2) - lambda * exp(-theta / 2) - theta / 2
```

The API is experimental and intentionally narrow:

```r
pc.penalty = list(scale = c(anchor, alpha))
pc.penalty = list(latent_scale = c(anchor, alpha))
```

The interpretation is `P(sigma_u > anchor) = alpha`; if `alpha` is omitted, it
defaults to `0.5`.

## Files touched

- `EBSmoothr/R/01_LGP.R`: LGP PC prior validation, R Laplace objective support,
  TMB data wiring, and returned diagnostics.
- `EBSmoothr/R/03_eb_smoother.R`: `family = "lgp"` wrapper support for
  `pc.penalty`, including learned-noise TMB fits.
- `EBSmoothr/src/EBSmoothr.cpp`: model_id 0 LGP TMB objective now optionally
  adds the latent-scale PC prior.
- `EBSmoothr/tests/testthat/test-lgp.R`: focused LGP PC prior tests.
- `EBSmoothr/man/ebnm_LGP_generator.Rd` and `EBSmoothr/man/eb_smoother.Rd`:
  draft documentation sync.

## Validation

Commands run:

```sh
Rscript --vanilla -e 'devtools::test("EBSmoothr", filter = "lgp", reporter = "summary")'
R CMD INSTALL -l /tmp/ebsmooth-lib EBSmoothr
Rscript --vanilla -e 'devtools::test("EBSmoothr", reporter = "summary")'
```

Both focused LGP tests and the full package test suite passed. The full suite
still emits INLA child-process retry/segfault messages in Matern tests, matching
previous behavior, but the R test run completed successfully.

## EBFPCA weather EBMF experiment

All EBMF checks used the draft package installed at `/tmp/ebsmooth-lib`, the
weather matrix from `EBFPCA/data/weather_data.rda`, and row-centered data:

```r
Yc <- Y - rowMeans(Y)
setup <- LGP_setup(t = seq_len(ncol(Yc)))
```

With `Kmax = 2` and `maxiter = 60`, plain LGP still hit the factor-2 iteration
cap:

| setting | factor 2 iterations | last ELBO diff | last LF max change |
| --- | ---: | ---: | ---: |
| plain | 60 | 7.17e-04 | 7.96e-08 |
| scale = c(1, 0.5) | 60 | 1.61e-03 | 7.88e-08 |
| scale = c(0.5, 0.5) | 60 | 1.61e-03 | 2.15e-07 |
| scale = c(0.25, 0.25) | 60 | 1.59e-03 | 4.69e-07 |

These broad anchors do not help because the fitted plain factor-2 latent scale
is much smaller: `theta ~= 11.7`, so `sigma_u ~= 0.0029`.

Near-scale anchors behave differently:

| setting | factor 2 iterations | last ELBO diff | last LF max change |
| --- | ---: | ---: | ---: |
| scale = c(0.005, 0.5) | 60 | 1.19e-03 | 2.59e-07 |
| scale = c(0.003, 0.5) | 60 | 9.63e-04 | 1.54e-07 |
| scale = c(0.002, 0.5) | 60 | 7.19e-04 | 1.72e-07 |
| scale = c(0.003, 0.1) | 8 | 1.73e-04 | 2.33e-05 |

The `scale = c(0.003, 0.1)` prior fixes the original factor-2 symptom in the
`Kmax = 2` reproduction, but it does not robustly solve the full greedy path.
With `Kmax = 5` and `maxiter = 30`, the same prior produced:

| factor | iterations | last ELBO diff | last LF max change | sigma_u |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 3 | 1.40e-05 | 8.14e-07 | 0.007552 |
| 2 | 8 | 1.73e-04 | 2.33e-05 | 0.002644 |
| 3 | 30 | 2.27e-03 | 2.10e-06 | 0.009099 |
| 4 | 30 | 1.60e-03 | 2.06e-06 | 0.006552 |
| 5 | 30 | 2.81e-03 | 4.37e-06 | 0.008112 |

Moving the anchor upward helps later factors but reintroduces factor-2
non-convergence. For example `scale = c(0.006, 0.1)` fixed factors 3-5 under
the same cap, but factor 2 hit 30 iterations.

## Draft conclusion

The implementation is numerically wired and passes tests, but the experiment
does not support adding this as a general convergence fix yet. A single global
latent-scale PC prior can locally distort the ridge when its anchor is matched
to a factor's fitted scale, but different EBMF factors land at different
`sigma_u` scales. This makes one common anchor a tuning-sensitive tradeoff
rather than a robust solution.

The observed behavior is still consistent with the scale-identifiability
diagnosis: `LF Max Chg` becomes tiny while ELBO continues moving. The PC prior
can change the ELBO slope along that direction, but it has to be calibrated to
the component scale to help.

## Posterior mean sanity check

A direct one-dimensional smoother check was added after the EBMF experiment. On
the same noisy curve, both `ebnm_LGP_generator()` and `eb_smoother()` produce
identical posterior means and fitted scales for a given PC prior. Stronger PC
priors shrink `sigma_u = exp(-theta / 2)` and reduce posterior mean roughness,
as expected.

For fixed `alpha = 0.1`, reducing the anchor from `5` to `0.5` to `0.03`
reduced latent scale and posterior roughness monotonically. For fixed
`anchor = 0.5`, reducing alpha from `0.9` to `0.5` to `0.02` also reduced latent
scale and posterior roughness monotonically.

This behavior is now covered in
`EBSmoothr/tests/testthat/test-lgp.R` by
`test_that("LGP PC prior smoothness responds to anchor and alpha", ...)`.

Additional validation:

```sh
Rscript --vanilla -e 'devtools::test("EBSmoothr", filter = "lgp", reporter = "summary")'
Rscript --vanilla -e 'devtools::test("EBSmoothr", reporter = "summary")'
```

Both passed. The full suite still prints INLA child-process retry/segfault
messages during Matern tests, matching existing behavior, but exits
successfully.

## Promotion to 0.2.6

The feature was promoted from draft exploration to an opt-in supported package
feature for EBSmoothr 0.2.6. The public API remains unchanged from the draft:
`pc.penalty = list(scale = c(anchor, alpha))` and
`pc.penalty = list(latent_scale = c(anchor, alpha))` are supported for L-GP
fits through `ebnm_LGP_generator()` and `eb_smoother(family = "lgp")`.

The documentation now describes the prior as supported scale regularization,
not as an experimental convergence fix. Matrix-factorization behavior remains a
motivating example rather than a package-level acceptance criterion.
