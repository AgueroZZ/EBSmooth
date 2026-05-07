# Matern Fisher-PQL Backend

## Summary

Implemented `backend = "fisher_pql"` for non-identity Matern fits and made it
the default `backend = "auto"` route for log and softplus Matern links after
inner-iteration sensitivity checks showed that three pseudo-Gaussian updates
retain most of the speedup while matching the Laplace reference closely on a
representative softplus `n = 1000` case. The first implementation incorrectly iterated full exact Gaussian
Matern fits and re-optimized range/sigma/noise/beta inside each PQL iteration.
The second implementation moved the one-step PQL solve into a TMB random-effect
Step A path, but runtime stayed close to identity-link TMB Laplace rather than
identity-link exact Gaussian. The current implementation uses the exact
pseudo-Gaussian Matern machinery for Fisher-PQL Step A.

## Changes

- Added Matern backend routing for `fisher_pql` with log/softplus-only
  validation.
- Added internal Fisher-PQL helpers for link derivatives, pseudo-response
  construction, sparse reference solves, exact Gaussian Step A fits, and
  PQL-mode final scoring.
- Added known-noise and learned-noise Fisher-PQL fit paths backed by the exact
  Matern solvers. Learned-noise Fisher-PQL passes a fixed `noise_scale` vector
  to the exact learned-noise solver.
- Removed the active TMB Fisher-PQL Step A route and the `model_id = 2` TMB
  template block.
- Kept `pql_inner_iter` internal and unexposed; the default is now three
  pseudo-Gaussian exact Matern updates.
- Changed the primary likelihood semantics to
  `fisher_laplace_at_fisher_pql_mode_<beta_mode>`. The reported
  `log_likelihood` is a Fisher/Laplace score evaluated at the PQL mode, not a
  true original-model Laplace marginal likelihood from re-optimizing the
  non-Gaussian latent field.
- Changed `backend = "auto"` for non-identity Matern fits to route to
  `fisher_pql`; identity-link auto routing remains exact Gaussian.
- Updated package docs and NEWS.

## Validation

- Parsed `EBSmoothr/R/02_Matern.R`, `EBSmoothr/R/03_eb_smoother.R`, and
  `EBSmoothr/tests/testthat/test-matern.R`.
- Regenerated roxygen documentation with `devtools::document("EBSmoothr")`.
- Ran `testthat::test_dir("EBSmoothr/tests/testthat", filter = "matern",
  stop_on_failure = TRUE)`: `424` passed, `0` failed, `0` skipped.
  Local INLA subprocesses emitted optimizer retry and crash messages in
  existing INLA/inlabru tests, but testthat completed successfully.
- Reinstalled the package with
  `R CMD INSTALL -l /Users/ziangzhang/Library/R/arm64/4.5/library EBSmoothr`;
  compilation completed successfully with Eigen unused-variable warnings.

## Runtime and Accuracy Notes

Exact Step A PQL now gives the intended runtime improvement on representative
learned-noise non-identity cases:

| n | link | Fisher-PQL | reference | speedup | identity exact | mean max/RMSE | eta max/RMSE | log range/sigma/noise error | objective drift / obs |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 300 | log | `1.311s` | `4.540s` | `3.46x` | `1.026s` | `0.157 / 0.042` | `0.089 / 0.030` | `0.045 / 0.036 / 0.025` | `-0.218` |
| 300 | softplus | `0.966s` | `3.442s` | `3.56x` | `1.026s` | `0.063 / 0.026` | `0.109 / 0.042` | `0.017 / 0.117 / 0.004` | `-0.075` |
| 1000 | log | `1.947s` | `23.993s` | `12.32x` | `2.011s` | `0.136 / 0.039` | `0.080 / 0.029` | `0.003 / 0.049 / 0.001` | `-0.122` |
| 1000 | softplus | `2.012s` | `22.832s` | `11.35x` | `2.011s` | `0.067 / 0.027` | `0.119 / 0.045` | `0.004 / 0.110 / 0.001` | `-0.064` |

The runtime gate passes strongly: Fisher-PQL is at least `3.4x` faster than
Laplace on these cases and is close to identity exact runtime for `n = 1000`.
Similarity is mixed. Hyperparameters, noise, eta-mode RMSE, and posterior RMSE
are close, but log-link posterior max error remains above `0.10`, and objective
drift per observation is materially larger than `0.02` for all four cases.
This backend should therefore remain opt-in and experimental rather than being
routed from `backend = "auto"`.

The existing exact learned-noise default path retained its fast path:

- `noise_scale = NULL`, `n = 80`: warmed timings `0.282, 0.281, 0.231s`;
  median `0.281s`.
- `noise_scale = rep(1, n)`, `n = 80`: warmed timings `0.237, 0.231,
  0.232s`; median `0.232s`.

Conclusion from the initial one-step validation: exact Step A Fisher-PQL is a
strong acceleration path for larger learned-noise non-identity Matern fits, but
one-step PQL was not uniformly close to Laplace under the objective and
max-error gates, especially for log-link learned-noise fits. This motivated the
inner-iteration sensitivity check below.

## Softplus Inner-Iteration Sensitivity

Follow-up experiment on `2026-05-07`: added an internal, non-public
`pql_max_iter` path for the Fisher-PQL helpers so we can test repeated exact
pseudo-Gaussian refits without changing the public API. Also stored
`posterior_eta` in Matern Laplace-style fit objects so latent-scale posterior
mean/sd diagnostics can be compared directly.

Script:

```sh
Rscript internal/simulations/benchmark_softplus_fisher_pql_inner_iter.R
```

Setup: `n = 1000`, two-dimensional softplus data, `max.edge = c(0.08, 0.24)`,
known-noise `ebnm_Matern_generator(..., backend = "laplace")` reference and
learned-noise `eb_smoother(..., s = NULL, backend = "laplace")` reference.
Fisher-PQL uses the exact pseudo-Gaussian Matern Step A with
`pql_max_iter = 1, 2, 3, 5, 10`.

| case | iter | PQL time | Laplace time | speedup | log-range diff | log-sigma diff | beta diff | noise SD diff | eta mean RMSE | eta sd RMSE | loglik diff / obs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| known-s ebnm | 1 | `4.228s` | `36.965s` | `8.74x` | `0.0151` | `-0.0735` | `0.0932` | NA | `0.1667` | `0.00687` | `-0.1428` |
| learned-s eb_smoother | 1 | `5.309s` | `38.679s` | `7.29x` | `0.0102` | `-0.0755` | `0.0929` | `-0.00067` | `0.1670` | `0.00685` | `-0.1587` |
| known-s ebnm | 2 | `6.302s` | `36.965s` | `5.87x` | `0.0014` | `-0.0100` | `0.0104` | NA | `0.0122` | `0.00087` | `-0.00081` |
| learned-s eb_smoother | 2 | `7.744s` | `38.679s` | `4.99x` | `-0.0017` | `-0.0115` | `0.0104` | `0.00014` | `0.0123` | `0.00081` | `-0.00090` |
| known-s ebnm | 3 | `8.391s` | `36.965s` | `4.41x` | `-0.0043` | `-0.0054` | `0.0063` | NA | `0.00028` | `0.00054` | `0.00009` |
| learned-s eb_smoother | 3 | `11.187s` | `38.679s` | `3.46x` | `-0.0047` | `-0.0056` | `0.0062` | `0.00002` | `0.00029` | `0.00055` | `0.00009` |
| known-s ebnm | 5 | `12.265s` | `36.965s` | `3.01x` | `-0.0045` | `-0.0051` | `0.0062` | NA | `0.00003` | `0.00054` | `0.00009` |
| learned-s eb_smoother | 5 | `14.961s` | `38.679s` | `2.59x` | `-0.0056` | `-0.0060` | `0.0061` | `0.00001` | `0.00003` | `0.00055` | `0.00010` |
| known-s ebnm | 10 | `22.483s` | `36.965s` | `1.64x` | `-0.0045` | `-0.0051` | `0.0062` | NA | `0.00004` | `0.00054` | `0.00009` |
| learned-s eb_smoother | 10 | `25.382s` | `38.679s` | `1.52x` | `-0.0056` | `-0.0060` | `0.0061` | `0.00001` | `0.00003` | `0.00055` | `0.00010` |

Result file:
`internal/simulations/results/softplus_fisher_pql_inner_iter_n1000.csv`.

Takeaway: one-step PQL is fastest but noticeably biased on this softplus case.
Two iterations remove most of the drift while retaining about `5x` speedup.
Three iterations are essentially indistinguishable from the Laplace reference
on latent-scale mean and objective diagnostics while retaining `3.5x-4.4x`
speedup. Iterations beyond three give negligible accuracy improvement and
mostly spend the speed advantage.

Based on this result, the default internal Fisher-PQL iteration count was
changed from `1` to `3`, and Matern `backend = "auto"` now routes log and
softplus links to `backend = "fisher_pql"`. Explicit `backend = "laplace"` and
`backend = "laplace_fisher"` remain available for direct reference fits.

Validation for this follow-up:

- Parsed `EBSmoothr/R/02_Matern.R` and
  `internal/simulations/benchmark_softplus_fisher_pql_inner_iter.R`.
- Ran the benchmark script above successfully.
- Ran an ad hoc non-INLA check covering `posterior_eta`, default public
  `backend = "fisher_pql"` dispatch, internal `pql_max_iter = 3`, and
  learned-noise Fisher-PQL output.
- After promoting Fisher-PQL to the non-identity Matern auto route, ran an
  ad hoc non-INLA check covering `ebnm_Matern_generator()` log/softplus auto
  routing, learned-noise `eb_smoother()` softplus auto routing, `posterior_eta`,
  and internal one-step `pql_max_iter = 1`.
- Regenerated roxygen documentation with `devtools::document("EBSmoothr")`.
- Reinstalled `EBSmoothr` version `0.2.3` with
  `R CMD INSTALL -l /Users/ziangzhang/Library/R/arm64/4.5/library EBSmoothr`.
- Built `EBSmoothr_0.2.3.tar.gz` with
  `R CMD build --no-build-vignettes EBSmoothr`.
- Ran `R CMD check --no-manual --no-build-vignettes --no-tests
  EBSmoothr_0.2.3.tar.gz`: completed with `Status: 1 NOTE`. The remaining
  NOTE is the existing `unlockBinding(default_name, ns)` safety note in
  `R/02_Matern.R`; newly introduced `setNames`, `.scale_w`, and `predict`
  notes were fixed.
- A full `testthat::test_dir(..., filter = "matern")` run was started but
  interrupted after the existing INLA/inlabru portion produced continuous
  `Hessian failed but no better mode found` diagnostics rather than reaching a
  testthat summary in reasonable time.
