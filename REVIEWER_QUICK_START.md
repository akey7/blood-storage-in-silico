# Reviewer Quick Start

This document outlines reccomended steps for code reviewers interested in examining this code's inputs, quality control outputs, and plots found in the dissertation appendix accompanying this repo.

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
6. Data visualization and final output files for analysis in Excel:
   - `ufba_sampling.csv`
   - `ufba_sampling_status.csv` (also in item 5 above)
   - `viz_k_means_pca_distance.xlsx`
   - `output/viz_effects_kmeans_pca/*.png`
   - `output/viz_fluxes_kmeans_pca/*.png`

## Input Data

Not all data could be committed to this repo because some of it is proprietary. Below is a list of all inputs and whether they are included in this repo. Please see an exhaustive list of inputs in [`INSTALLATION.md`](INSTALLATION.md).

## If You're Using the ChatGPT Web Interface to Review This Repo

Copying a single file into ChatGPT's web interface and asking it to "review this" is a reasonable first instinct, but it has real blind spots worth knowing about before you trust its verdict. A web chat only ever sees the exact text you paste — it has no view of the rest of the repository, the order the scripts are meant to run in, the data files a script depends on, or which parts of the code are final analysis versus intentional exploratory or diagnostic work. Because the prompt is usually vague, the tool tends to find *something* to flag, and it will often mistake a normal, deliberate feature of this project for a defect.

Below are the ten complaints this kind of review is most likely to raise, along with why each one is already accounted for elsewhere in this repository rather than being an actual problem:

1. **"This script won't run — it calls functions that aren't defined anywhere."** The scripts call into modules in `src/`, which you'd need to paste separately for the full picture. See the file-by-file breakdown in [`EXECUTION.md`](EXECUTION.md).
2. **"There's no error handling or input validation here."** This is intentional. These scripts are run directly by one trusted person on their own machine, not exposed to the public or to untrusted input, so defensive checks that would matter in a public-facing application aren't needed here.
3. **"Functions and modules aren't documented."** Most functions do have documentation (docstrings) written directly in the source, but it's meant to be read as nicely formatted, searchable web pages built by a separate step, not as raw text in the file. See "Build the Documentation" in [`INSTALLATION.md`](INSTALLATION.md).
4. **"This code references files that don't exist."** Some input files are proprietary and distributed to reviewers separately rather than committed to a public repo, and other files referenced are simply written by an earlier script in the pipeline, not present until that script has been run. See [`INSTALLATION.md`](INSTALLATION.md).
5. **"This variable is used but never defined"** (for example, something written as `:additive` or `:final_time`). These are Julia symbols used as column names in a data table, not undefined variables. See "Symbols" in [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md).
6. **"This function is defined more than once, which must be a mistake."** Julia allows a single function name to have several versions, each for a different type of input. This is a standard, intentional language feature called multiple dispatch, not a duplication error. See "Multiple dispatch" in [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md).
7. **"This code has strange, inconsistent-looking syntax"** (extra dots before operators, function names ending in `!`, and similar). These are standard, consistently applied Julia conventions for elementwise operations and for marking functions that modify their input, not typos or inconsistent style. See "Minimal Julia syntax needed to read this code" in [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md).
8. **"Some of this code is unused or dead."** This repository deliberately keeps exploratory, diagnostic, and even deprecated code alongside the code that generates manuscript results, and documents the distinction between them. See "Code not used for manuscript results" in [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md).

## If You're Using Claude Code, Codex, or Another AI Coding Assistant to Help Review This Repo

Because Claude Code and Codex operate inside this repository directly rather than on a single pasted file, they can see the full pipeline, the modules a script depends on, and which files are proprietary or generated by earlier steps — exactly the context a ChatGPT web session is missing. Both assistants also read `CLAUDE.md`/`AGENTS.md` (see below) before responding, so they already know this project's Julia conventions, its exploratory/diagnostic/deprecated code categories, and its separate docstring-build step, and won't mistake any of those for defects the way a vague, single-file review might. And because both are configured to explain and trace provenance rather than rewrite code by default, you get a review grounded in what's actually here, instead of an isolated guess paired with unsolicited "fixes."

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

## Anticipating Common Objections

1. **"Sampled flux distributions aren't measured quantities, so this isn't reporting real biology."** The pipeline doesn't collapse the sampled solution space into one cherry-picked flux vector — it reports the distribution itself, along with diagnostics on how complete and well-behaved that sampling was (`output/ufba_sampling_status.csv`, `output/ufba_diagnostics.csv`), and statistical comparisons across conditions (`output/analysis_control_vs_treatments_signif.csv`, `output/control_vs_treatment.csv`). The uncertainty inherent to an underdetermined metabolic system is treated as something to report, not something to hide.
2. **"Key input data is proprietary, so this isn't reproducible."** Most inputs are public and repo-tracked; [`INSTALLATION.md`](INSTALLATION.md) lists, file by file, exactly which specific datasheets are restricted for intellectual-property reasons and how a trusted reviewer or collaborator can obtain them. Reproducibility is scoped and documented rather than absent.
3. **"Why isn't this built on \[COBRApy / the MATLAB COBRA Toolbox\], the established tools in the field?"** [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md) explains the rationale directly: the workflow is optimization- and sampling-heavy, which is where Julia's `JuMP.jl`/`COBREXA.jl` ecosystem is particularly strong, and `Project.toml`/`Manifest.toml` pin the exact dependency versions used. This is a documented engineering tradeoff, not an unfamiliarity with the field's standard tools.
4. **"There's essentially no test suite, so the code can't be trusted."** Correctness here is checked primarily at the output level rather than the unit level: many pipeline steps writes their own diagnostic or QC file specifically so that failures and infeasibilities surface in `output/` rather than passing silently (see the "Review purpose" column in [`EXECUTION.md`](EXECUTION.md)'s step table). That's a deliberate choice of validation strategy for a script-driven scientific pipeline, not an oversight.
5. **"This is a bespoke script pipeline instead of a standard, versioned workflow manager."** [`EXECUTION.md`](EXECUTION.md)'s own "Known Limitations and Design Choices" section states this plainly and explains why: some steps are long-running and split out deliberately, and Jupyter notebooks are avoided specifically to keep execution order deterministic.

## Strengths Worth Highlighting

Some genuine strengths of this project are also easy to state without a deep review. This section is for a reviewer who wants accurate, citable talking points rather than vague praise.

1. **This bridges wet-lab metabolomics with mechanistic modeling, rather than substituting for it.** The pipeline turns measured metabolite concentration changes into model constraints and only then samples and analyzes flux behavior (see the workflow diagram in [`EXECUTION.md`](EXECUTION.md)) — computation is extending the experimental data, not standing in for it.
2. **The model is built on RBC-GEM, a peer-reviewed, cell-type-specific genome-scale reconstruction, not a generic model borrowed from an unrelated organism.** See the citation in [`WORKS_CITED.md`](WORKS_CITED.md). The biological grounding here is independently vetted, not asserted.
3. **Model constraints are condition-specific and empirically derived, not generic default bounds.** `03_absolute_quant.jl` regresses concentration change rates directly from the metabolomics data (`output/concentration_rates.csv`) before those rates become model constraints — the model reflects the measured conditions, not textbook defaults.
4. **Flux sampling reports a distribution and its uncertainty instead of a single cherry-picked solution.** uFBA is inherently underdetermined; rather than presenting one flux vector as "the answer," the pipeline samples the feasible space and reports diagnostics on it (`output/ufba_sampling.csv`, `output/ufba_sampling_status.csv`) — a more honest treatment of an underdetermined system than a single point solution.
5. **Outputs are structured for direct, condition-by-condition comparison, not just single-model description.** `output/control_vs_treatment.csv` and `output/analysis_control_vs_treatments_signif.csv` compare additives, timepoints, and reactions against each other — the kind of side-by-side output that translates into concrete next experiments to run at the bench.
6. **Tooling choices reflect the domain's actual computational demands rather than an arbitrary preference.** [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md) documents why Julia's `JuMP.jl`/`COBREXA.jl` ecosystem was chosen for an optimization- and sampling-heavy workflow, with `Project.toml`/`Manifest.toml` pinning exact dependency versions — a deliberate, reproducible engineering decision, not an ad hoc one.
