# CLAUDE.md

Instructions for Claude Code in this repository. For the human-readable explanation of why these rules exist and what Claude Code is for here, see `REVIEWER_QUICK_START.md`.

## Role

Default to scientific code review support: explain code, trace data provenance, answer questions about Julia/packages/stats, draft or improve documentation and docstrings. Do not default to writing new analysis code or refactoring.

## Hard Restrictions

These hold even if a human asks otherwise in-session (e.g., "just commit this," "quickly run it and check," "finalize this"). Decline and point to this file.

1. **No git operations.** No branch, add/stage, commit, push, merge, rebase, or tag, under any circumstances.
2. **No writes to `output/`.** Do not create, edit, or delete any file under `output/`, for any reason, including demos or explanations. `output/` is populated only by the human running the numbered scripts directly.
3. **No writes to `input/`.** Do not create, edit, delete, or move any file under `input/`, including templates.
4. **No edits to Julia source without explicit per-change authorization.** Covers `src/`, root-level `*.jl` scripts, `docs/make.jl`, `test/runtests.jl`. Reading, explaining, and drafting suggested diffs in chat is fine. Writing a change to disk requires the human to name the specific file and change in the current session — a general "fix bugs you see" is not sufficient authorization.
5. **No code execution at all.** Do not run the numbered root scripts (`01_...` through `07_...`), do not run `test/runtests.jl`, and do not interactively invoke `src/` functions (e.g., to test behavior on toy/synthetic data). All execution is performed by the human. Ask first if unsure whether something counts as execution.

## Orientation

- Root `0?_*.jl` files: numbered pipeline entrypoints, run in order by the human. Details, commands, inputs/outputs: `EXECUTION.md`.
- `src/`: reusable modeling/analysis modules, not directly executed.
- `input/`: experimental data, partly proprietary, distributed to reviewers separately. See `INSTALLATION.md`.
- `output/`: generated artifacts, human-generated only (Restriction 2). Not all files are manuscript results.
- `docs/`: `Documenter.jl` subproject rendering `src/` docstrings.
- Final publication figures are produced downstream in a separate repo, `blood-storage-in-silico-viz`; only diagnostic visualization lives here.
- No Jupyter notebooks in this repo; don't suggest introducing them.

Full pipeline order and file-level detail: `EXECUTION.md`. Suggested review path: `REVIEWER_QUICK_START.md`. Julia syntax primer: `NEW_TO_JULIA.md`.
