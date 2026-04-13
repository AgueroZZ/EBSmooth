# Refresh Matern Internal Simulations

## What changed

- Reframed the 1D parameter-recovery study as an increasing-domain experiment with fixed sampling density.
- Standardized the 1D sample sizes at 80, 160, 320, and 640.
- Kept the 2D validation simple but stronger than before by using larger grids and a higher observation noise level.
- Refreshed the internal report so the CSV summaries, figures, and markdown/HTML outputs come from the same script version.

## Why

The earlier internal validation already showed good surface recovery, but it mixed older and newer result files and the 2D example was visually too easy. This refresh is meant to provide a cleaner convergence check in 1D and a more informative noisy-surface comparison in 2D.

## Key findings

- In the 1D increasing-domain study, the mean fitted range moved from 0.1729 at `n = 80` to 0.1954 at `n = 640`, with truth 0.2000.
- In the same 1D study, the mean fitted sigma moved from 0.2527 at `n = 80` to 0.2949 at `n = 640`, with truth 0.3000.
- In the 2D recovery study, the larger grid (`36 x 28`, `n = 1008`) gave mean fitted range 0.2164 and mean fitted sigma 0.3545 for truths 0.2200 and 0.3500.
- The noisier 2D example (`30 x 24`, `noise sd = 0.15`) produced a visibly rougher observed surface while the posterior mean still tracked the latent field well (`RMSE = 0.1078`, `corr = 0.9731`).
- The profile marginal-likelihood surface for the 2D example looked plausible: the optimized point and the true point lay on the same high-likelihood region, with the true point 2.0501 log-likelihood units below the coarse-grid maximum.
- The exact marginal-likelihood implementation still matched a dense Gaussian calculation to machine precision, with maximum absolute difference `5.68434e-14`.
