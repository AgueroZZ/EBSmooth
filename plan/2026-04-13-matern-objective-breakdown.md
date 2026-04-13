## Objective

Decouple the expensive exact-objective cross-check from the main INLA PC-prior
Matern fit, and provide a post-fit function that can recompute the exact
Gaussian objective components on demand.

## Plan

1. Turn the automatic exact cross-check in the non-fixed PC-prior Matern path
   into an explicit option, with the default set to off.
2. Store enough Matern context in fitted objects so that exact objectives can
   be recomputed after fitting.
3. Add an exported helper that reports the exact fixed-beta, profiled-beta,
   and integrated-flat objectives, along with the PC-prior contribution when
   relevant.
4. Update tests and internal documentation to use the new post-fit workflow.
