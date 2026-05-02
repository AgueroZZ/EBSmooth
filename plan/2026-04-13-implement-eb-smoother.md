# Implement `eb_smoother()` for Known and Learned Noise

## Summary
- Add a new public API `eb_smoother()` to support both known observation standard errors and learned common noise SDs.
- Keep `ebnm_LGP_generator()` and `ebnm_Matern_generator()` unchanged as `ebnm`-compatible low-level interfaces.
- Implement the first version for both `matern` and `lgp`, with `inla_pc` support only on the Matern path.

## Key Changes
- Add `eb_smoother()` as a unified user-facing entrypoint.
- Add learned-noise exact and INLA-PC Matern fitting paths.
- Add learned-noise LGP fitting by extending the TMB objective with a scalar noise parameter.
- Return a new fit class `eb_smoother_fit` with familiar fields plus `fitted_noise_sd` and `noise_mode`.
- Update the EBMF experiment helpers to call `eb_smoother()` for posterior-mean post-smoothing.

## Tests
- Regression checks for known-`s` fits against current Matern and LGP implementations.
- Learned-noise checks for exact Matern, INLA-PC Matern, and LGP.
- API validation for unsupported family/backend combinations and invalid `pc.penalty$noise` usage.

## Assumptions
- V1 learns only one common scalar noise SD when `s = NULL`.
- V1 does not expose learned-noise fits as `ebnm` objects.
- Matern remains identity-link only; LGP keeps `identity` and `log`.
