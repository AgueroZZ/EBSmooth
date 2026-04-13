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
- Each individual fit is capped at `120` seconds.

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
completed | 27
timed_out | 8

Plot: `/Users/ziangzhang/Desktop/EBSmooth/internal/simulations/results/runtime_scaling_1d_methods.png`

## Timing Table

n | lgp_b0 | status_lgp_b0 | lgp_blt0 | status_lgp_blt0 | matern_b0_exact | status_matern_b0_exact | matern_blt0_exact | status_matern_blt0_exact | matern_pc_b0_inla | status_matern_pc_b0_inla | matern_pc_blt0_exact | status_matern_pc_blt0_exact | matern_pc_blt0_clinear | status_matern_pc_blt0_clinear
--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---
100.0000 | 0.2770 | completed | 0.4160 | completed | 0.9680 | completed | 0.8620 | completed | 2.9870 | completed | 0.7570 | completed | 3.0200 | completed
300.0000 | 0.6050 | completed | 0.5770 | completed | 0.8280 | completed | 2.4300 | completed | 3.6010 | completed | 3.0010 | completed | 3.7560 | completed
1000.0000 | 7.0430 | completed | 6.6230 | completed | 0.7190 | completed | 12.2130 | completed | 3.4150 | completed | 8.1420 | completed | 3.3480 | completed
3000.0000 | 120.0000 | timed_out | 120.0000 | timed_out | 0.8090 | completed | 120.0000 | timed_out | 4.3010 | completed | 120.0000 | timed_out | 3.3300 | completed
10000.0000 | 120.0000 | timed_out | 120.0000 | timed_out | 2.7300 | completed | 120.0000 | timed_out | 17.4480 | completed | 120.0000 | timed_out | 4.4720 | completed
