# Reviewer Quick Start

This document outlines reccomended steps for code reviewers interested in examining this code's inputs, quality control outputs, and outputs for downstream data visualization.

## Suggested Review Path

For a scientific review, inspect the code in this order:

1. [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md): Enough Julia to understand the broad strokes of this repo if you are coming from a Python or R background.
2. [`INSTALLATION.md`](INSTALLATION.md): How to install the dependencies in this repo.
1. [`EXECUTION.md`](EXECUTION.md) to understand the workflow order.
2. `03_absolute_quant.jl` and `src/AbsoluteQuant.jl` to understand how metabolomics data are converted into concentration change rates.
3. `04_ufba_sampler.jl` and the following modules to understand how the uFBA sampling is performed:
   - `src/FbaModelBuilder.jl`
   - `src/MetaboliteBounds.jl`
   - `src/PruningOptimizations.jl`
   - `src/UfbaSampler.jl`
4. Downstream flux summaries and statistical comparisons: 
   - `05_ufba_sampler_analysis.jl` and `src/UfbaSamplerAnalysis.jl`
   - `06_ufba_sampler_analysis_2.jl` and `src/UfbaSamplerAnalysis2.jl`
5. Diagnostic outputs in `output/`, especially:
   - `ufba_sampling_status.csv`
   - `ufba_pruning_overview.csv`
   - `ufba_fba_breaks.csv`
   - `ufba_prune_breaks.csv`
   - `measurements_and_sinks_report.csv`
   - `measurements_and_sinks_report_by_model.csv`
6. Outputs the drive data visualization in `blood-storage-in-silico-viz` repo:
   - `ufba_sampling.csv`
   - `ufba_sampling_status.csv` (also in item 5 above)
   - `viz_k_means_pca_distance.xlsx`

## Input Data

Not all data could be committed to this repo because some of it is proprietary. Below is a list of all inputs and whether they are included in this repo. Please see an exhaustive list of inputs in [`INSTALLATION.md`](INSTALLATION.md).

## If You're Using Claude Code, Codex, or Another AI Coding Assistant to Help Review This Repo

This repository includes a `CLAUDE.md` file that configures Claude Code's behavior here, and an `AGENTS.md` file that gives the same configuration to Codex and other agents that look for that filename instead. `AGENTS.md` defers to `CLAUDE.md` as the canonical source of truth and only adds notes for applying the rules in a Codex-style workflow, so the two files enforce identical restrictions regardless of which tool you're running. This section explains, in plain language, what that configuration means and why it exists — whether you're the one running the assistant, or you're a reviewer reading this repository and want to know what an AI assistant was and wasn't allowed to do.

**The intended role of Claude Code, Codex, or any other assistant in this repository is scientific code review support — not code generation.** Concretely, that means it's set up to default to:

- Explaining what a script or function does, and why, in plain language.
- Helping you trace data provenance: which input file feeds which script, which script produces which output file, and how a given manuscript figure or table traces back to raw data.
- Answering questions about Julia syntax, package usage, or statistical/modeling choices (see `NEW_TO_JULIA.md` if you're coming from Python or R).
- Drafting or improving standalone documentation, such as README-style `.md` files, when asked. Docstrings and comments live inside Julia source, so instead of drafting or editing those directly, it will explain which ones are missing, unclear, or out of date (see Restriction 4 below).

It is explicitly configured **not** to write new analysis code, refactor existing modules, run any code, or touch git, by default. Peer review of this repository is about verifying what was actually run to produce the manuscript's results — not about producing a different or "improved" analysis. If you see the assistant declining a request along these lines, that's expected behavior, not a malfunction.

### Restrictions Claude Code and Codex follow in this repo

1. **No git operations.** It will not create branches, stage files, commit, push, merge, or tag, even if asked to "save" or "finalize" work. All git operations are performed manually by the repository maintainer.
2. **No writes to `output/`.** It will not create, edit, or delete anything in `output/`, including during a demo or explanation. That folder is populated only by the human running the numbered scripts directly.
3. **No writes to `input/`.** It will not create, edit, delete, or move any file in `input/`, including template files. Some input files are proprietary and distributed to reviewers through separate, controlled channels (see `INSTALLATION.md`).
4. **No edits to Julia source, and no suggested diffs, under any circumstances.** This covers `src/` modules (including docstrings and comments), the numbered root scripts, `docs/make.jl`, and `test/runtests.jl`. This cannot be lifted in-session — no per-change sign-off, explicit naming of the file and change, or any other request from a human user makes writing to these files permitted, and none of it permits producing a diff, patch, or full replacement snippet either, even one explicitly framed as "just for the human to paste in." It may still read and explain this code, and describe in prose what a fix or a missing docstring would involve.
5. **No code execution, period.** It will not run the numbered pipeline scripts, run `test/runtests.jl`, or interactively invoke functions from `src/` — including just to demonstrate behavior on toy data. All execution, including test runs, is performed by the human. This is because the numbered scripts must run in a fixed order, some are long-running and resource-intensive, and they depend on proprietary inputs and a fully configured `output/` directory (see `INSTALLATION.md` and `EXECUTION.md`).

If you want to see exactly how these rules are phrased for the assistant itself, they're in `CLAUDE.md` at the repo root for Claude Code, and `AGENTS.md` at the repo root for Codex and similar tools.
