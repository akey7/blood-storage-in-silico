# Blood Storage *In Silico*

This project studies refrigerated blood storage with unsteady flux balance analysis (uFBA).

**Too long, didn't read:** Reviewers should see the [REVIEWER_QUICK_START.md](REVIEWER_QUICK_START.md) guide for information on code structure, setup, suggested review path, and AI governance in tailored to your use case.

## FAQ

**What does this repository actually do?**
It models refrigerated blood storage using unsteady flux balance analysis (uFBA): metabolomics measurements are converted into metabolic constraints, fed into a red blood cell genome-scale model, and sampled to produce flux distributions that are then analyzed statistically. See the step-by-step workflow in [EXECUTION.md](EXECUTION.md).

**Why is this written in Julia instead of Python or R?**
The core workflow is optimization- and sampling-heavy — model construction, feasibility checks, and constrained flux sampling — which is where Julia's optimization ecosystem (`JuMP.jl`, `COBREXA.jl`) is particularly strong. [NEW_TO_JULIA.md](NEW_TO_JULIA.md) explains the rationale in more detail and includes a syntax primer if you're coming from Python or R.

**I don't know Julia. Can I still review this code?**
Yes. [NEW_TO_JULIA.md](NEW_TO_JULIA.md) is written for exactly this situation — it covers the syntax patterns and packages you'll encounter without requiring you to become a Julia programmer.

**Where should I start if I want to review the science?**
[REVIEWER_QUICK_START.md](REVIEWER_QUICK_START.md) lays out a recommended reading order through the scripts, modules, and diagnostic outputs, ending with the files that feed the downstream visualization repo.

**Can I run the pipeline myself?**
Yes, if you have Julia and the required input files. [INSTALLATION.md](INSTALLATION.md) covers installing dependencies and setting up `input/`/`output/`, and [EXECUTION.md](EXECUTION.md) gives the exact commands for each numbered script, in order.

**Why can't I access all of the input data?**
Some metabolomics datasheets are proprietary and distributed to trusted reviewers and collaborators separately rather than committed to this public repo. [INSTALLATION.md](INSTALLATION.md) lists exactly which files are public versus restricted, and how to obtain the restricted ones.

**Is there a required order for running the scripts?**
Yes — the root-level `0?_*.jl` scripts are numbered entrypoints meant to run in that sequence, since later steps consume the outputs of earlier ones. The full dependency chain, including which script writes which files, is diagrammed in [EXECUTION.md](EXECUTION.md).

**Where do the final published figures come from?**
Diagnostic and intermediate visualizations are generated in this repo, but the publication-quality figures are produced downstream in a separate repository, `blood-storage-in-silico-viz`.

**Does this repo use Jupyter notebooks?**
No, intentionally. Notebooks are avoided in favor of the numbered scripts to keep execution order deterministic and reproducible; see [EXECUTION.md](EXECUTION.md) for the reasoning.

**Are AI coding assistants like Claude Code or Codex allowed to work in this repo?**
Only in a constrained, review-support role — explaining code, tracing data provenance, and drafting documentation, but not writing new analysis code, running scripts, or touching git by default. These constraints are configured in [CLAUDE.md](CLAUDE.md) and [AGENTS.md](AGENTS.md) and explained in plain language in [REVIEWER_QUICK_START.md](REVIEWER_QUICK_START.md).

**What papers underpin the modeling approach?**
See [WORKS_CITED.md](WORKS_CITED.md) for the key references, including the RBC-GEM genome-scale model and the metabolomics platform used to generate the underlying data.

## How to Get Started

Start with [REVIEWER_QUICK_START.md](REVIEWER_QUICK_START.md): it lays out a recommended order for reviewing this code, and explains how [CLAUDE.md](CLAUDE.md) and [AGENTS.md](AGENTS.md) constrain Claude Code and Codex, respectively, to the same scientific-review-support role in this repo — read it first if you plan to use either assistant while reviewing.

Once you're ready to set up the repo, [INSTALLATION.md](INSTALLATION.md) covers installing Julia and its dependencies, obtaining the (partly proprietary) input files, and building the `src/` docstring documentation. [EXECUTION.md](EXECUTION.md) then walks through the numbered pipeline scripts in order, including the known limitations and design choices behind how they're split up.

If you're new to Julia, [NEW_TO_JULIA.md](NEW_TO_JULIA.md) gives enough context on syntax and packages to follow the scientific logic without becoming a Julia programmer. For sources underpinning the modeling choices, see [WORKS_CITED.md](WORKS_CITED.md).

## Quick Hints: Mental Model of this Repository

This repository is organized as a script-driven Julia analysis pipeline.

- Root-level `*.jl` files are command-line entrypoints.
- Files in `src/` contain the reusable modeling, data-processing, optimization, sampling, analysis, and diagnostic visualization logic.
- Root-level scripts mostly connect `src/` modules to the filesystem, command-line options, and the ordered workflow.
- The primary scientific workflow is documented in `EXECUTION.md`.
- Installation, proprietary input setup, and output directory setup are documented in `INSTALLATION.md`.
- Publication-quality final figures are produced in a separate Python/R repository, `blood-storage-in-silico-viz`.
