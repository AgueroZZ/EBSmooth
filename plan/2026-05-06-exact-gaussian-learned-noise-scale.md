# Exact Gaussian Learned-Noise Scale Extension

## Goal

Extend the internal exact-Gaussian learned-noise Matern solver so future
Fisher-PQL pseudo-response fits can estimate one scalar noise SD while using a
fixed per-observation scale:

```r
s_i = noise_sd * noise_scale_i
```

For Fisher-PQL, `noise_scale_i = 1 / g_i`, where `g_i` is fixed inside one
Gaussian solve and updated only between PQL iterations.

## Implementation Plan

- Add an internal helper that maps scalar `noise_sd` and optional
  `noise_scale` to the effective observation SD vector.
- Keep `noise_scale = NULL` on the existing fast path, using the current
  `rep(noise_sd, n)` construction.
- Thread `noise_scale` through all exact learned-noise objective modes:
  empirical-Bayes beta, fixed beta, flat beta prior, and proper beta prior.
- Preserve PC-prior semantics on scalar `noise_sd`, not on effective
  per-observation SDs.
- Store the scalar `fitted_noise_sd` and effective `fitted_s` in the exact
  learned-noise fit object.

## Validation Plan

- Verify `noise_scale = NULL` remains the default exact learned-noise path.
- Test `noise_scale = rep(1, n)` against `noise_scale = NULL`.
- Test heterogeneous `noise_scale` across all beta modes.
- Compare heteroskedastic learned-noise fits against exact known-noise objective
  evaluations at `s = fitted_noise_sd * noise_scale`.
- Run a local timing smoke check before and after the change on the default
  `noise_scale = NULL` path.
