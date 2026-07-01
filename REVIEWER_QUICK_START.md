# Reviewer Quick Start

This document outlines reccomended steps for code reviewers interested in examining this code's inputs, quality control outputs, and outputs for downstream data visualization.

## Suggested Review Path

For a scientific review, inspect the code in this order:

1. [`NEW_TO_JULIA.md`](NEW_TO_JULIA.md): Enough Julia to understand the broad strokes of this repo if you are coming from a Python or R background.
2. [`INSTALLATION.md`](INSTALLATION.md): How to install the dependencies in this repo.
1. [`EXECUTION.md`](EXECUTION.md) to understand the workflow order.
2. `absolute_quant.jl` and `src/AbsoluteQuant.jl` to understand how metabolomics data are converted into concentration change rates.
3. `ufba_sampler.jl` and the following modules to understand how the uFBA sampling is performed:
   - `src/FbaModelBuilder.jl`
   - `src/MetaboliteBounds.jl`
   - `src/PruningOptimizations.jl`
   - `src/UfbaSampler.jl`
4. Downstream flux summaries and statistical comparisons: 
   - `ufba_sampler_analysis.jl` and `src/UfbaSamplerAnalysis.jl`
   - `ufba_sampler_analysis_2.jl` and `src/UfbaSamplerAnalysis2.jl`
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
