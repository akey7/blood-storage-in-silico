# Commands to Execute Scripts

## Known Limitations and Design Choices

- Root-level scripts are used as workflow entrypoints instead of a single workflow manager.
- Some workflow steps are split into separate scripts because they are long-running.
- Proprietary input files are not committed to the public repository.
- `output/` contains generated artifacts and diagnostic files; not every file in `output/` is a manuscript result.
- Jupyter notebooks are intentionally avoided to preserve deterministic script execution order.

## Multithreaded Execution

Note that in both macOS commands, the `JULIA_NUM_THREADS` environment variable in the execution command sets the number of threads avaiable to Julia.

The commands to launch the scripts on Windows are similar, but the JULIA_NUM_THREADS environment variable is not specified in the command line. Rather, the environment variable is configured in settings as described in [INSTALLATION.md](INSTALLATION.md).

Either way, customize the `JULIA_NUM_THREADS` environment variable to match your environment.

## Setting Up `input/` and `output/` Folders

Before you execute these scripts, the `input/` and `output/` folder must be configured as described in [INSTALLATION.md](INSTALLATION.md).

## uFBA Workflow: Order of Script Execution

### Table of Steps

Conceptually, the modeling process involves these steps:

```text
experimental metabolomics data
        ↓
condition-specific metabolite constraints
        ↓
red blood cell metabolic model
        ↓
feasibility checks and constraint construction
        ↓
optimization steps
        ↓
flux sampling
        ↓
summary tables of sampled flux distributions
        ↓
statistical analysis and visualization
```

This table offers a summary of the order of script execution. Full details of script execution are below.

Note: **Steps 3-6 are the core workflow, while steps 1, 2, and 7 are optional.**

| Step | Entrypoint                     | Main module(s)                                                                                           | Main inputs                                    | Main outputs                             | Review purpose                             |
| ---: | ------------------------------ | -------------------------------------------------------------------------------------------------------- | ---------------------------------------------- | ---------------------------------------- | ------------------------------------------ |
|    1 | `01_plot_metabolite_timelines.jl` | `src/MetaboliteTimelines.jl`                                                                             | relative metabolomics inputs                   | timeline plots, correlations             | QC of metabolite trajectories              |
|    2 | `02_raw_relative_intensities.jl`  | `src/RawRelativeIntensities.jl`                                                                          | relative quant data                            | PCA loadings and PCA plots               | exploratory QC                             |
|    3 | `03_absolute_quant.jl`            | `src/AbsoluteQuant.jl`                                                                                   | relative + absolute quant data                 | concentration rates                      | converts metabolomics to model constraints |
|    4 | `04_ufba_sampler.jl`              | `src/UfbaSampler.jl`, `src/FbaModelBuilder.jl`, `src/PruningOptimizations.jl`, `src/MetaboliteBounds.jl` | concentration rates, RBC-GEM, opt-in/out files | uFBA models, sampled fluxes, diagnostics | core modeling/sampling step                |
|    5 | `05_ufba_sampler_analysis.jl`     | `src/UfbaSamplerAnalysis.jl`                                                                             | sampled fluxes                                 | median fluxes, diagnostics, comparisons  | primary analysis outputs                   |
|    6 | `06_ufba_sampler_analysis_2.jl`          | `src/UfbaSamplerAnalysis2.jl`                                                                                  | sampled fluxes, median fluxes                  | diagnostic plots, PCA/k-means spreadsheet  | Final plots in dissertation appendix             |
|    7 | `07_model_graph.jl`               | `src/ModelGraph.jl`                                                                                      | uFBA models, DFS plan                          | graph traversal outputs                  | graph-based inspection                     |

### Conceptual Graphical Map

```text
[Raw/Relative Quant MS Data] ──> (1) 01_plot_metabolite_timelines.jl (Diagnostic plots of raw relative intensity data)
                                 (2) 02_raw_relative_intensities.jl (PCAs raw relative intensity data)
                                      │
[Absolute Quant Datasheets]  ───> (3) 03_absolute_quant.jl ──> [Rates & Concentrations]
                                                              │
[RBC-GEM Metabolic Model]    ─────────────────────────────────┴─> (4) 04_ufba_sampler.jl (Long-running sampling)
                                                                       │
                                      ┌────────────────────────────────┘
                                      ▼
                                 (5) 05_ufba_sampler_analysis.jl (Long-running analysis)
                                 (6) ufba_sampler_analysis_2.jl (Faster analysis)
                                 (7) 07_model_graph.jl (Graph/DFS Traversal)
```

### Map of Significant Data Outputs

This is a non-exhaustive diagram of files that are generated on `output/` by the various steps of the workflow. Even though it is non-exhaustive, it does list the files most important for analysis and reporting results from the modeling run.

```text
input metabolomics files
        │
        ▼
03_absolute_quant.jl
        │
        └── output/concentration_rates.csv
        │
        ▼
04_ufba_sampler.jl
        │
        ├── output/ufba_sampling.csv
        ├── output/ufba_sampling_status.csv
        ├── output/ufba_pruning_overview.csv
        └── output/ufba_models/
        │
        ▼
05_ufba_sampler_analysis.jl
        │
        ├── output/ufba_median_fluxes.csv
        ├── output/analysis_control_vs_treatments_signif.csv
        ├── output/measurements_and_sinks_report.csv
        └── output/flux_vector_data_matrices/
ufba_sampler_analysis_2.jl
        │
        ├── output/ufba_sampling_complete_additives.csv
        └── output/viz_k_means_pca_distance.xlsx
        └── output/viz_effects_kmeans_pca/
        └── output/viz_fluxes_kmeans_pca/
```

## *Optional*: Execute the First 2 Pipeline Steps to Explore Relative Quantification Data

All commands are issued from the root of the repo.

### (1) `01_plot_metabolite_timelines.jl`: Plot Timelines of Relative Metabolite Intensities

Uses `src/MetaboliteTimelines.jl` to generate plots of timerseries of metabolite timelines.

On macOS and Linux

```
JULIA_NUM_THREADS=7 julia --project=. 01_plot_metabolite_timelines.jl
```

on Windows (assuming your `JULIA_NUM_THREADS` environment variable is set)

```
julia --project=. 01_plot_metabolite_timelines.jl
```

Output will be saved to `output/normalized_abundance_correlations.csv` and `output/plots`

### (2) `02_raw_relative_intensities.jl`: PCA Plots of Relative Quant Data

Uses `src/RawRelativeIntensities.jl` to make PCA plots reducing relative metabolite abundances down to fewer features.

Execution is multithreaded, so the number of threads should be specified.

On macOS and Linux:

```
JULIA_NUM_THREADS=7 julia --project=. 02_raw_relative_intensities.jl
```

On Windows, assuming `JULIA_NUM_THREADS` has been set in settings:

```
julia --project=. 02_raw_relative_intensities.jl
```

Outputs:
1. PCA loadings for all additives to `output/relative_pca_loadings.csv`.
2. PCA plot DataFrames as `.csv` files to `output/pca_plot_dfs`.
3. 2D PCA plots of single additives and pairs of additives to `output/pca_plots`.

## *Required*: Execute the Core Workflow (Steps 3 to 6)

### (3) `03_absolute_quant.jl`: Approximate Absolute Quantifications and Regress Concentration Change Rates

Uses `src/AbsoluteQuant.jl` to perform the following tasks:

1. Loads data relative (`Data Sheet 1.CSV`) and absolute quant (`Absolute Quant Data Sheet.xlsx` and `Absolute Quant Extracellular Datasheet.xlsx`) data files.

2. Writes quality checks to `output/qc_fold_changes.csv` and `output/qc_fold_change_zeros.csv`.

3. Map and proportinoate (see documentation built in [INSTALLATION.md](INSTALLATION.md)) combined metaboilite names to single metabolite ids from the RBC-GEM. Approximate absolute concentrations for the relative quant data using accompanying absolute quant data. Writes result to `output/relative_absolute_quant.csv`.

4. Plots the approximated absolute concentrations over time to plots in the `output/relative_absolute_plots/[cleaned metabolite name].png`

5. Perform c-means clustering on the metabolite trajectories. Create plots of c-means clusters and an accompanying elbow plot for each additive to `output/relative_absolute_c_means/`. Write cluster memberships to `output/c_means_primary_clusters.csv`.

6. Performs the regression to determine the rates of metabolite concentration changes and writes the result to `output/concentration_rates.csv`. Plots regression data and stores the plots in `output/regression_plots/`

This script uses multiple threads to calculate all the regression quickly, so it relies on the `JULIA_NUM_THREADS` variable.

On macOS and Linux, execute the following (customize the number of threads to your machine):

```
JULIA_NUM_THREADS=7 julia --project=. 03_absolute_quant.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. 03_absolute_quant.jl
```

### (4) `04_ufba_sampler.jl`: Run uFBA Sampling Jobs

In addition to multithreading, the uFBA sampling module uses concurrent worker processes to fully utilize the hardware executing the script. There is an optimal point to set the number of workers: if there are too few workers, the job will take a needlessly long time to execute. With too many workers, the script takes too long to launch.

**Default to 4 workers with `-p 4` below unless you have a good reason to do otherwise.** Greater numbers of workers will lengthen script startup time.

**At the end of sampling, job status will be displayed**. Not every model will execute. That isn't a script error; rather, it is a shortcoming of input data. Subsequent data processing will step over missing data. More detail is explained in the paper.

**uFBA sampling may take some time**: The uFBA sampling is a computationally intensive step and the duration of the task will depend on the number of cores you have devoted to the task with the command line options and the number of threads specified with the `JULIA_NUM_THREADS` environment variable.

The command line arguments to the Julia environment and script are the following:

1. `-p`: How many workers Julia will lauch for the sampling task.

2. `--nchains`: The number of chains of samples for each uFBA model. More chains means more sampling and longer execution.

3. `--nmodels`: Number of uFBA models to analyze (-1 for all possible models).

4. `--prune-method`: The pruning method to use. Cane be either `case1` or `case3`. See Bordbar (2016) or `PruningOptimizations.jl` for more information.

On a macOS or Linux machine with 14 cores, an example command to set the number of workers and threads on the same line would be (while executing all models with 5 chains and case3 pruning):

```
JULIA_NUM_THREADS=7 julia --project=. -p 4 04_ufba_sampler.jl --nchains 5 --nmodels -1 --prune-method case3
```

On a Windows machine, an example to work with your previously set `JULIA_NUM_THREADS` environment variable would be (again while executing all models with 5 chains and case3 pruning):

```
julia --project=. -p 4 04_ufba_sampler.jl --nchains 5 --nmodels -1 --prune-method case3
```

Which would sample all models with 5 chains, run all models, and use 4 concurrent workers.

Customize workers, threads, number of chains, and number of models your use case. For quick runs, set the number of models and chains to be small numbers.

In addition to input files from prior steps, there is an input file of note
1. `input/metabolite_measurement_opt_outs.csv`: A file of metabolite ids of absolute quant measurements (without the leading `M_`) to ignore when building all models. Used for diagnostic purposes for failing models.
2. `input/sink_opt_ins.csv`: Metabolite ids that should have sinks. This list overrides pruning decisions. If this file is missing and you need a template for it, see `input/sink_opt_ins_template.csv`, which is tracked in source control.

Outputs the following files:
1. `output/ufba_sampling_status.csv`: That statuses of each uFBA sampling job (fail or ok)
2. `output/ufba_sampling.csv`: The samplings of the fluxes. Used by next step.
4. `output/fba_model_metabolites.csv`: Metabolite ids of the FBA models created for the uFBA runs.
5. `output/ufba_blocked_reactions.csv`: Reaction ids of blocked reactions and their corresponding strings for each model.
6. `output/ufba_optimized_sinks.csv`: Reaction ids of sinks sampled, their corresponding metabolites and directions, and median fluxes.
7. `output/debug_case1.lp` (if configured in the code): Diagnostic output from `optimize_case_1()` to assist in debugging Case 1 optimization runs.
8. `output/debug_case3.lp` (if configured in the code): Diagnostic output from `optimize_case_3()` to assist in debugging Case 1 optimization runs.
9. `output/ufba_prune_breaks.csv`: Constraints broken in pruning attempts across all uFBA models.
10. `output/ufba_fba_breaks.csv`: Constraints broken in simple FBA attempts executed before the uFBA runs.
11. `output/ufba_pruning_overview.csv`: Zero and non-zero sinks found in the pruning process. Helpful to see what decisions the pruning algorithm made.

### (5) `05_ufba_sampler_analysis.jl`: Analyze the results of the uFBA Runs

Outputs `.csv` and `.xlsx` analyses of the uFBA results. These files can be used by themselves and they also feed into the next step of making visualizations.

Runs code in the `src/UfbaSamplerAnalysis.jl`.

Uses the following input files:

1. Reads the uFBA sampling results file at `output/ufba_sampling.csv`.

Outputs the following files:

1. Diagnoses the output of the models sampled by uFBA to help find potential problems and writes the diagnostics in `output/ufba_diagnostics.csv`.
2. Writes net fluxes of each pair of sinks to `output/net_sink_fluxes.csv`.
3. Writes data matrices of median fluxes to `output/flux_vector_data_matrices`. One matrix contains all additives. The rest of the matrices exclude one matrix at a time.
4. Writes a report of all metabolites in each model and whether those metabolites are measured or have sinks to `output/measurements_and_sinks_report.csv`.
5. Writes an aggregated report for each model detailing the total numbers of metabolites, measurements, and sinks to `output/measurements_and_sinks_report_by_model.csv`.
6. `output/control_vs_treatment.csv`: Potentially interesting additives/times/reactions for further investigation. See the documentation for the function `compare_flux_distributions()` in `UfbaSamplerAnalysisAndViz.jl` for more information.
7. Writes a bunch of `.csv` files for use by an R script to plot correlation heatmaps and perform hierarchical clustering. These files are written to `output/correlation_matrices_1/`.
8. Writes an Excel workbook that links reactions to metabolites and counts the number of measures metabolties per reaction, reaction subsystem, and reaction category. Filename is `output/reactions_metabolites_measurements.xlsx`
9. Writes an Excel workbook of comparing treatments with respect to reactions and time points. Filename is `output/reaction_treatment_comparison.xlsx`

On macOS and Linux, set the `JULIA_NUM_THREADS` environment variable and execute like this:

```
JULIA_NUM_THREADS=7 julia --project=. 05_ufba_sampler_analysis.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. 05_ufba_sampler_analysis.jl
```

**Note**: Because the most time-consuming part of the analysis is multithreaded, dots are printed instead of a green status indicator for thread safety status reporting.

### (6) `06_ufba_sampler_analysis_2.jl`: Visualize uFBA analysis results as plots

Creates visualizations (histograms, KDE plots, PCA scatter plots) of the uFBA analysis results.

Uses the following input files:

1. Reads the uFBA sampling results file at `output/ufba_sampling.csv`.
2. Reads the reactions ids to strings YAML file at `output/rxn_ids_to_strings.yml`.
3. Reads the median reaction fluxes at `output/ufba_median_fluxes.csv`.

Outputs the following files:

1. Writes histograms of sampling results (one plot per reaction) to `output/uFBA_histograms_v2/`.
2. Writes kernel density estimation of sampling results (one plot per reaction) to `output/uFBA_densities/`.
3. PCA and k-means analysis of Cohen's effects between control and treatments and median fluxes, along with distances between control and treatments to sheets in `output/viz_k_means_pca_distance.xlsx`.

On macOS and Linux, set the `JULIA_NUM_THREADS` environment variable and execute like this:

```
JULIA_NUM_THREADS=7 julia --project=. 06_ufba_sampler_analysis_2.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. 06_ufba_sampler_analysis_2.jl
```

### (7) `07_model_graph.jl`: Analyze the uFBA models as graphs

Analyzes the uFBA models as graphs.

Runs code in `src/ModelGraph.jl`.

The "DFS plan" file in `input/gem_dfs/dfs_plan.csv` needs the following columns:

1. `metabolite_id`: Metabolite id to start a DFS traversal at.
2. `max_depth`: Maximum number of hops to traverse from the starting vertex.

A default `input/dfs_plan.csv` is provided in the repo as an example.

Inputs
1. `input/dfs_plan.csv`, which are the metabolites to start traversal of the graph from.
2. uFBA sampling results file at `output/ufba_sampling.csv`.
3. Reaction id to reaction string mapping file at `output/rxn_ids_to_strings.yml`.
4. uFBA model SBML files in `output/ufba_models`.

Outputs
1. Writes `output/gem_dfs/visited_metabolites.csv` (which specifies the metabolites traversed on DFS traversals)
2. Writes `output/gem_dfs/visited_reaction.csv` (which specifies the reactions traversed on DFS traversals).

Displays progress bars to show progress as it works through the data.

This script only uses a single thread, so execution On macOS and Linux or Window is simple:

```
julia --project=. 07_model_graph.jl
```

## *Optional*: Tests

There is a small test to ensure testing works:

On macOS and Linux:

```
julia --project=. test/runtests.jl
```

On Windows

```
julia --project=. .\test\runtests.jl
```
