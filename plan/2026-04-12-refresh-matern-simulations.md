# Refresh Matern Internal Simulations

## Goal

Refresh the internal Matern validation so that it better addresses two statistical questions:

1. In 1D, does the empirical-Bayes hyperparameter estimate behave plausibly under an increasing-domain design?
2. In 2D, does the fitted posterior mean visibly separate the latent smooth surface from noisier observed data?

## Planned Changes

- Keep the 1D setup simple but switch the interpretation to an increasing-domain regime with roughly fixed sampling density.
- Use 1D sample sizes 80, 160, 320, and 640 to provide a clearer convergence trend without making the run too expensive.
- Use a larger and noisier 2D setup than the earlier low-noise example, while keeping the mesh and profile-likelihood grid modest.
- Re-render the figures, CSV summaries, and markdown/HTML report from the same script so the outputs are internally consistent.

## Deliverables

- Updated `internal/simulations/run_matern_validation.R`
- Refreshed files under `internal/simulations/results/`
- A short project log describing the new settings and the refreshed conclusions
