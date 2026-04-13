# Plan: Release Readiness and Publish Pass

Date: 2026-04-13

## Goal

Do a final release-readiness sweep for the current `EBSmoothr` state, tighten
documentation and package housekeeping, verify the package checks cleanly, and
publish the resulting changes to the remote repository.

## Steps

1. Inspect the current repository state, package metadata, and documentation to
   identify publish-facing gaps.
2. Add a repository-level README that explains the package scope, installation,
   and quick-start usage.
3. Clean package metadata and build rules so generated artifacts and hidden
   files are not shipped in source builds.
4. Fix the R-package MIT license stub while preserving the full license text in
   a companion markdown file.
5. Remove local junk files and compiled artifacts from the working tree.
6. Run package verification checks on the built source package.
7. Stage, commit, and push the final release-readiness changes.
