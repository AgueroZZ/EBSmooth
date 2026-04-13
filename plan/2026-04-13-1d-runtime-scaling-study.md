## Objective

Run a lightweight 1D scalability study across the main LGP and Matern
algorithmic variants, covering `betaprec = 0` and `betaprec < 0`, and save a
scatterplot of runtime versus sample size.

## Plan

1. Define a simple 1D benchmark design with increasing sample size and a fixed
   noise level.
2. Include the following variants:
   - LGP with `betaprec = 0`;
   - LGP with `betaprec < 0`;
   - Matern exact with integrated-flat beta (`betaprec = 0`);
   - Matern exact with profiled beta (`betaprec < 0`);
   - Matern + PC prior with `betaprec = 0`;
   - Matern + PC prior with `betaprec < 0`, exact/profile implementation;
   - Matern + PC prior with `betaprec < 0`, `clinear` implementation.
3. Benchmark one fit per method per sample size and save the raw timing table.
4. Generate a log-log runtime plot and a short internal markdown summary.
5. Record the study in `log/` for future collaborators.
