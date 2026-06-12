# Draft Plan: LGP Latent-Scale PC Prior

## Goal

Add an experimental LGP-only `pc.penalty` option that places an exponential
PC prior on the LGP latent standard deviation `sigma_u = exp(-theta / 2)`.
Use it to test whether a scale-aware prior mitigates flashier's slow ELBO
drift along the loading/factor scale ridge in the EB-FPCA weather example.

## Design

- Support `pc.penalty = list(scale = c(anchor, alpha))`.
- Also accept `pc.penalty = list(latent_scale = c(anchor, alpha))` as a more
  explicit alias.
- Interpret the anchor as `P(sigma_u > anchor) = alpha`, where
  `sigma_u = exp(-theta / 2)`.
- Add the log prior to the optimized LGP marginal objective, not as a post-hoc
  log-likelihood adjustment.
- Keep the feature off by default and preserve existing behavior when
  `pc.penalty = NULL`.

## Validation

- Add focused tests for parsing, log-prior contribution, invalid arguments,
  `ebnm_LGP_generator()`, and `eb_smoother(family = "lgp")`.
- Rebuild the draft package in a temporary library.
- Repeat the weather EBMF convergence case with and without the LGP PC prior.
- Compare iteration behavior, reconstruction change, and factor scaling drift.
