---
title: "Inspect INLA Step-A and Step-B Objective Quantities"
---

# Inspect INLA Step-A and Step-B Objective Quantities

Generated on 2026-04-12 16:29:04 CDT

## Setup

- Fixed single 2D dataset on a 20 x 16 grid (n = 320).
- Truth: range = 0.30, sigma = 0.35, beta0 = 0.20, noise sd = 0.22, alpha = 2.
- Three legacy INLA scenarios are compared by changing the PC-prior tail probabilities only.

## Step A vs Step B

scenario | range_prob | sigma_prob | fitted_range | fitted_sigma | fitted_beta | stepA_mlik_integration | stepA_mlik_gaussian | stepA_log_posterior_mode | stepA_max_log_posterior | stepA_joint_log_posterior | exact_loglik_at_stepA_mode | stepB_mlik_integration | stepB_mlik_gaussian
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
balanced | 0.500000 | 0.500000 | 0.303481 | 0.362997 | 0.195837 | -97.002754 | -94.877195 | -92.516545 | -92.516545 | -96.715072 | -89.052058 | -89.052058 | -89.052058
smooth_lowvar | 0.100000 | 0.100000 | 0.310747 | 0.364053 | 0.195600 | -97.871755 | -95.746196 | -93.360016 | -93.360017 | -97.584073 | -89.082946 | -89.082946 | -89.082946
rough_highvar | 0.900000 | 0.900000 | 0.300609 | 0.362631 | 0.195931 | -99.569059 | -97.443500 | -95.090844 | -95.090844 | -99.281377 | -89.048650 | -89.048650 | -89.048650

## Step B With Fixed Theta/Beta Across Different Priors

scenario | range_prob | sigma_prob | stepB_fixed_theta_mlik_integration | stepB_fixed_theta_mlik_gaussian | stepB_fixed_theta_log_posterior_mode | stepB_fixed_theta_max_log_posterior
--- | --- | --- | --- | --- | --- | ---
balanced | 0.500000000000 | 0.500000000000 | -89.052058053207 | -89.052058053207 | 0.000000000000 | -89.052058053207
smooth_lowvar | 0.100000000000 | 0.100000000000 | -89.052058053207 | -89.052058053207 | 0.000000000000 | -89.052058053207
rough_highvar | 0.900000000000 | 0.900000000000 | -89.052058053207 | -89.052058053207 | 0.000000000000 | -89.052058053207

## Interpretation

- `stepB_mlik_*` matches the exact Gaussian marginal likelihood evaluated at the fitted hyperparameters.
- `stepB_mlik_*` and `stepB_fixed_theta_*` are invariant to the PC-prior choice once theta is fixed.
- The Step A quantities that change with the prior are `stepA_mlik_*`, `stepA_log_posterior_mode`, `stepA_max_log_posterior`, and `stepA_joint_log_posterior`.
- Empirically, `stepA_log_posterior_mode` and `stepA_max_log_posterior` behave like penalized optimization targets, whereas Step B does not retain a prior-dependent criterion once theta is fixed.
