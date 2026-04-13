## Update

Added a lightweight 2D comparison study between exact Matern smoothing and an
estimated iid normal prior under two truths:

- a true 2D Matern field;
- a true iid Gaussian field.

## Completed Changes

- Added `internal/simulations/study_2d_matern_vs_iid_comparison.R`.
- Added two settings:
  - `smooth_truth`: latent surface simulated from the exact Matern model;
  - `iid_truth`: latent surface simulated independently across locations.
- For each replicate, compared:
  - exact Matern EB fit;
  - `ebnm::ebnm_normal(mode = "estimate", scale = "estimate")`.
- Saved per-replicate draws and separate summaries for:
  - model selection;
  - smooth-truth hyperparameter recovery;
  - iid-truth fitted range behavior.
- Added summary figures for selection probabilities and hyperparameter behavior.

## Validation

- `Rscript internal/simulations/study_2d_matern_vs_iid_comparison.R --smoke`
- `Rscript internal/simulations/study_2d_matern_vs_iid_comparison.R`

## Notes

- The study is intentionally lightweight: two 2D grids and a small number of
  replicates per cell.
- The goal is to provide a readable empirical comparison, not a formal
  asymptotic statement.

## Lightweight Run

- Default run configuration:
  - grids: `18 x 14` and `24 x 18`
  - reps per cell: `6`
  - true Matern hyperparameters: `range = 0.22`, `sigma = 0.35`, `beta0 = 0.2`
  - iid latent scale: `tau0 = 0.35`
  - noise sd: `0.15`
- Main output artifacts:
  - `internal/simulations/results/2d_matern_vs_iid_comparison_draws.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_selection_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_smooth_truth_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_iid_truth_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_summary.md`
- Headline findings:
  - Under the true 2D Matern setting, the marginal likelihood selected the
    smooth Matern model in all `12/12` replicates.
  - In the same smooth setting, fitted hyperparameters stayed close to the
    truth and improved on the larger grid.
  - Under the iid setting, the marginal likelihood split `6/12` for iid and
    `6/12` for Matern in this lightweight run.
  - When Matern won under iid truth, the fitted range remained very small
    relative to the mesh edge, consistent with a near-boundary fit rather than
    genuine smooth-structure recovery.

## PC-Prior Extension

- Extended `internal/simulations/study_2d_matern_vs_iid_comparison.R` with an
  optional PC-prior mode controlled by:
  - `--use_pc_prior=true`
  - `--pc_range_anchor`
  - `--pc_range_alpha`
  - `--pc_sigma_anchor`
  - `--pc_sigma_alpha`
- In PC-prior mode, the script passes `pc.penalty` to
  `ebnm_Matern_generator(...)` and writes results under the
  `2d_matern_vs_iid_comparison_pc_*` prefix by default.

## PC-Prior Validation

- `Rscript internal/simulations/study_2d_matern_vs_iid_comparison.R --use_pc_prior=true --smoke`
- `Rscript internal/simulations/study_2d_matern_vs_iid_comparison.R --use_pc_prior=true`

## PC-Prior Lightweight Run

- PC-prior configuration:
  - range penalty: `P(range < 0.22) = 0.1`
  - sigma penalty: `P(sigma > 0.35) = 0.5`
- Main output artifacts:
  - `internal/simulations/results/2d_matern_vs_iid_comparison_pc_draws.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_pc_selection_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_pc_smooth_truth_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_pc_iid_truth_summary.csv`
  - `internal/simulations/results/2d_matern_vs_iid_comparison_pc_summary.md`
- Headline findings:
  - Under the true 2D Matern setting, the PC-prior run still selected the
    smooth Matern model in all `12/12` replicates.
  - Under the iid setting, the PC-prior run selected the iid model in all
    `12/12` replicates, compared with `6/12` in the no-prior baseline.
  - The paired replicate comparison showed that all six iid-truth replicates
    that previously selected Matern switched to iid after turning on the
    PC prior.
  - Under smooth truth, fitted ranges and sigmas moved slightly toward the
    true values and the selection margin remained strongly in favor of the
    smooth model.
  - Under iid truth, the penalized Matern fits no longer collapsed to
    near-zero ranges; instead, the PC prior pulled the fitted range upward
    while making the penalized objective decisively worse than the iid model.
