# Plan: Clean Commit Documentation Synchronization

Date: 2026-05-03

## Goal

Prepare the repository for a clean documentation-only commit after auditing the
current package state.

## Planned Changes

- Keep untracked scratch artifacts out of the next package commit by using local
  git excludes rather than deleting local experiment files.
- Update internal vignette sources so they describe the current public API:
  `family = "constant"`, Fisher Laplace log-link defaults, exact identity-link
  PC-prior auto fits, and explicit INLA compatibility aliases.
- Re-render internal vignette HTML outputs from the updated sources.
- Re-run source parse, Rd checks, testthat, and git status before reporting the
  final commit scope.

## Commit Scope

- Include only synchronized internal vignette sources, rendered internal HTML
  outputs, and this plan/log record.
- Do not include `.vscode/`, `attempts_in_EBMF/`, `tryingSpatial.R`, compiled
  objects, local session files, or generated attempt figures.
