# Blood Storage *In Silico*

This project studies refrigerated blood storage with unsteady flux balance analysis (uFBA).

## Documentation

Start with [REVIEWER_QUICK_START.md](REVIEWER_QUICK_START.md): it lays out a recommended order for reviewing this code, and explains how [CLAUDE.md](CLAUDE.md) constrains Claude Code (or another AI coding assistant) to a scientific-review-support role in this repo — read it first if you plan to use an assistant while reviewing.

Once you're ready to set up the repo, [INSTALLATION.md](INSTALLATION.md) covers installing Julia and its dependencies, obtaining the (partly proprietary) input files, and building the `src/` docstring documentation. [EXECUTION.md](EXECUTION.md) then walks through the numbered pipeline scripts in order, including the known limitations and design choices behind how they're split up.

If you're new to Julia, [NEW_TO_JULIA.md](NEW_TO_JULIA.md) gives enough context on syntax and packages to follow the scientific logic without becoming a Julia programmer. For sources underpinning the modeling choices, see [WORKS_CITED.md](WORKS_CITED.md).

## Reviewer Mental Model

This repository is organized as a script-driven Julia analysis pipeline.

- Root-level `*.jl` files are command-line entrypoints.
- Files in `src/` contain the reusable modeling, data-processing, optimization, sampling, analysis, and diagnostic visualization logic.
- Root-level scripts mostly connect `src/` modules to the filesystem, command-line options, and the ordered workflow.
- The primary scientific workflow is documented in `EXECUTION.md`.
- Installation, proprietary input setup, and output directory setup are documented in `INSTALLATION.md`.
- Publication-quality final figures are produced in a separate Python/R repository, `blood-storage-in-silico-viz`.
