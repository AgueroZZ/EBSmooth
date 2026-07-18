# Standalone EBSmoothr repository migration

## Objective

Create `AgueroZZ/EBSmoothr` as the canonical public R package repository while
retaining `AgueroZZ/EBSmooth` as the research and writing workspace.

## Steps

1. Audit the current `EBSmoothr/` subtree, its commit history, and any nested
   dirty worktrees before changing repository topology.
2. Extract package history with `git subtree split --prefix=EBSmoothr`.
3. Create the standalone GitHub repository and push the extracted history to
   its `main` branch.
4. Add package-facing metadata, a user README with a copy-paste `pak`
   installation block, and continuous integration appropriate for an R
   package repository.
5. Run package tests, build, and `R CMD check` from the standalone repository.
6. Tag and publish `v0.2.6` in the standalone repository.
7. Replace `EBSmooth/EBSmoothr` with a submodule only after all local content
   is safely published or backed up, then update the research-repository
   README to point users to the canonical package repository.
8. Verify both remote repositories and record the final commit, tag, release,
   and submodule SHAs in `log/`.

## Safety constraints

- Preserve unrelated untracked experiments, writing, output, and plans.
- Do not stage the mixed working tree wholesale.
- Do not delete or overwrite nested worktree changes.
- Use explicit paths for every staged commit.
