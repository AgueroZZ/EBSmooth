# 1D Runtime Scaling Study

This internal benchmark studies how runtime scales with sample size `n`
for the main 1D algorithmic variants currently available or internally
prototyped in this repository.

## Design

- Sample sizes: `100, 300, 1000, 3000, 10000`.
- Data are generated on a 1D increasing domain with approximately 100
  observation points per unit length.
- Noise level is fixed at `0.22`.
- The runtime includes the full end-to-end fit for each method, including
  setup/generator construction.
- For Matern methods, the mesh edge is fixed at `0.2` to let latent-field
  complexity grow with the domain length.
- For L-GP, the knot count is chosen from the domain length with a minimum
  of 30 and a maximum of 500.
- Each individual fit is capped at `60` seconds.

## Methods in the plot

- `LGP`, `betaprec = 0`
- `LGP`, `betaprec < 0`
- `Matern`, `betaprec = 0`, exact integrated-flat beta objective
- `Matern`, `betaprec < 0`, exact profiled beta objective
- `Matern + PC prior`, `betaprec = 0`, INLA Step A with beta integrated out
- `Matern + PC prior`, `betaprec < 0`, exact profiled objective + PC prior
- `Matern + PC prior`, `betaprec < 0`, INLA `clinear` prototype

## Status counts

status | count
--- | ---
completed | 29
error | 4
timed_out | 2

Plot: `/Users/ziangzhang/Desktop/EBSmooth/internal/simulations/results/runtime_scaling_1d_methods.png`

## Timing Table

n | lgp_b0 | status_lgp_b0 | lgp_blt0 | status_lgp_blt0 | matern_b0_exact | status_matern_b0_exact | matern_blt0_exact | status_matern_blt0_exact | matern_pc_b0_inla | status_matern_pc_b0_inla | matern_pc_blt0_exact | status_matern_pc_blt0_exact | matern_pc_blt0_clinear | status_matern_pc_blt0_clinear
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
100.0000 | 0.1400 | completed | 0.1660 | completed | 1.4150 | completed | 0.4990 | completed | 2.5730 | completed | 0.8400 | completed | 2.7190 | completed
300.0000 | 0.2400 | completed | 0.2130 | completed | 0.9420 | completed | 0.6090 | completed | 3.0030 | completed | 2.6290 | completed | 2.8690 | completed
1000.0000 | 2.2330 | completed | 2.1230 | completed | 0.7980 | completed | 0.9350 | completed | 3.0100 | completed | 21.7650 | completed | 3.0190 | completed
3000.0000 | 60.0000 | timed_out | 60.0000 | timed_out | 0.7200 | completed | 2.9190 | completed | 3.0330 | completed | 60.1360 | error | 2.9150 | completed
10000.0000 | 163.1080 | error | 160.6370 | error | 2.1160 | completed | 18.0900 | completed | 5.1280 | completed | 60.0630 | error | 4.1290 | completed

## Public API Benchmark at `n = 3000`

method_id | elapsed_seconds | status
--- | --- | ---
lgp_known | 60.0000 | timed_out
lgp_learned | 60.0000 | timed_out
matern_exact_known | 2.5530 | completed
matern_exact_learned | 2.4850 | completed
matern_pc_known | 3.4300 | completed
matern_pc_learned | 3.8640 | completed

## Repeated Matern Smoothing Benchmark at `n = 3000`

method_id | n_components | elapsed_seconds | status
--- | --- | --- | ---
rebuild_each_call | 4.0000 | 12.2250 | completed
reuse_matern_setup | 4.0000 | 11.2730 | completed
