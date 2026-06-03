# Commands to Execute Scripts

## Multithreaded Execution

Note that in both macOS commands, the `JULIA_NUM_THREADS` environment variable in the execution command sets the number of threads avaiable to Julia.

The commands to launch the scripts on Windows are similar, but the JULIA_NUM_THREADS environment variable is not specified in the command line. Rather, the environment variable is configured in settings as described in [INSTALLATION.md](INSTALLATION.md).

Either way, customize the `JULIA_NUM_THREADS` environment variable to match your environment.

## Setting Up `input/` and `output/` Folders

Before you execute these scripts, the `input/` and `output/` folder must be configured as described in [INSTALLATION.md](INSTALLATION.md).

## uFBA Workflow - Order of Script Execution

The scripts take input and write output files. Some scripts rely on output files previously written by other scripts. The order of script execution presented here maintains the order of reliance of the scripts on each other if such order is important. Such an arrangment of scripts might not be ideal, but it works for now!

## Manually Executing Scripts

All commands are issued from the root of the repo.

### (1) `plot_metabolite_timelines.jl`: Plot Timelines of Relative Metabolite Intensities

Uses `src/MetaboliteTimelines.jl` to generate plots of timerseries of metabolite timelines.

on macOS

```
JULIA_NUM_THREADS=7 julia --project=. plot_metabolite_timelines.jl
```

on Windows (assuming your `JULIA_NUM_THREADS` environment variable is set)

```
julia --project=. plot_metabolite_timelines.jl
```

Output will be saved to `output/normalized_abundance_correlations.csv` and `output/plots`

### (2) `raw_relative_intensities.jl`: PCA Plots of Relative Quant Data

Uses `src/RawRelativeIntensities.jl` to make PCA plots reducing relative metabolite abundances down to fewer features.

Execution is multithreaded, so the number of threads should be specified.

On macOS:

```
JULIA_NUM_THREADS=7 julia --project=. raw_relative_intensities.jl
```

On Windows, assuming `JULIA_NUM_THREADS` has been set in settings:

```
julia --project=. raw_relative_intensities.jl
```

Outputs:
1. PCA loadings for all additives to `output/relative_pca_loadings.csv`.
2. PCA plot DataFrames as `.csv` files to `output/pca_plot_dfs`.
3. 2D PCA plots of single additives and pairs of additives to `output/pca_plots`.

### (3) `absolute_quant.jl`: Approximate Absolute Quantifications and Regress Concentration Change Rates

Uses `src/AbsoluteQuant.jl` to perform the following tasks:

1. Loads data relative (`Data Sheet 1.CSV`) and absolute quant (`Absolute Quant Data Sheet.xlsx` and `Absolute Quant Extracellular Datasheet.xlsx`) data files.

2. Writes quality checks to `output/qc_fold_changes.csv` and `output/qc_fold_change_zeros.csv`.

3. Map and proportinoate (see documentation built in [INSTALLATION.md](INSTALLATION.md)) combined metaboilite names to single metabolite ids from the RBC-GEM. Approximate absolute concentrations for the relative quant data using accompanying absolute quant data. Writes result to `output/relative_absolute_quant.csv`.

4. Plots the approximated absolute concentrations over time to plots in the `output/relative_absolute_plots/[cleaned metabolite name].png`

5. Perform c-means clustering on the metabolite trajectories. Create plots of c-means clusters and an accompanying elbow plot for each additive to `output/relative_absolute_c_means/`. Write cluster memberships to `output/c_means_primary_clusters.csv`.

6. Performs the regression to determine the rates of metabolite concentration changes and writes the result to `output/concentration_rates.csv`. Plots regression data and stores the plots in `output/regression_plots/`

This script uses multiple threads to calculate all the regression quickly, so it relies on the `JULIA_NUM_THREADS` variable.

On macOS, executethe following (customize the number of threads to your machine):

```
JULIA_NUM_THREADS=7 julia --project=. absolute_quant.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. absolute_quant.jl
```

### (4) `ufba_sampler.jl`: Run uFBA Sampling Jobs

In addition to multithreading, the uFBA sampling module uses concurrent worker processes to fully utilize the hardware executing the script. There is an optimal point to set the number of workers: if there are too few workers, the job will take a needlessly long time to execute. With too many workers, the script takes too long to launch.

The command line arguments to the Julia environment and script are the following:

1. `-p`: How many workers Julia will lauch for the sampling task.

2. `--nchains`: The number of chains of samples for each uFBA model. More chains means more sampling and longer execution.

3. `--nmodels`: Number of uFBA models to analyze (-1 for all possible models).

4. `--prune-method`: The pruning method to use. Cane be either `case1` or `case3`. See Bordbar (2016) or `PruningOptimizations.jl` for more information.

On a macOS or Linux machine with 14 cores, an example command to set the number of workers and threads on the same line would be (while executing all models with 5 chains and case3 pruning):

```
JULIA_NUM_THREADS=7 julia --project=. -p 4 ufba_sampler.jl --nchains 5 --nmodels -1 --prune-method case3
```

On a Windows machine, an example to work with your previously set `JULIA_NUM_THREADS` environment variable would be (again while executing all models with 5 chains and case3 pruning):

```
julia --project=. -p 4 ufba_sampler.jl --nchains 5 --nmodels -1 --prune-method case3
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

### (5) `ufba_sampler_analysis.jl`: Analyze the results of the uFBA Runs

Outputs `.csv` and `.xlsx` analyses of the uFBA results. These files can be used by themselves and they also feed into the next step of making visualizations.

Runs code in the `src/UfbaSamplerAnalysis.jl`. Shows nifty status bars to indicate progress.

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

On macOS, set the `JULIA_NUM_THREADS` environment variable and execute like this:

```
JULIA_NUM_THREADS=7 julia --project=. ufba_sampler_analysis.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. ufba_sampler_analysis.jl
```

### (6) `ufba_sampler_viz.jl`: Visualize uFBA analysis results as plots

Creates visualizations (histograms and KDE plots) of the uFBA analysis results.

Uses the following input files:

1. Reads the uFBA sampling results file at `output/ufba_sampling.csv`.
2. Reads the reactions ids to strings YAML file at `output/rxn_ids_to_strings.xml`

Outputs the following files:

1. Writes histograms of sampling results (one plot per reaction) to `output/uFBA_histograms_v2/`.
8. Writes kernel density estimation of sampling results (one plot per reaction) to `output/uFBA_densities/`.

On macOS, set the `JULIA_NUM_THREADS` environment variable and execute like this:

```
JULIA_NUM_THREADS=7 julia --project=. ufba_sampler_viz.jl
```

On Windows, ensure that `JULIA_NUM_THREADS` is set and execute:

```
julia --project=. ufba_sampler_viz.jl
```

### (7) `model_graph.jl`: Analyze the uFBA models as graphs

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

This script only uses a single thread, so execution on macOS or Window is simple:

```
julia --project=. model_graph.jl
```

## API Runner

The `api_runner.jl` script unifies many aspects of the manual workflow above for execution by the Python API interface.

## Other Scripts

There are other scripts that you can run in this project. They are outside of the main uFBA workflow, and are thus optional. They are documented here for completeness.

### `test_case_1.jl` and `test_case_3.jl`: Test Case 1 and Case 3 sink pruning

I built this script to test Case 1 and Case 3 sink pruning code and to serve as an example for more involved workflows in `UfbaSampler.jl`.

Case 1 on macOS:

```
JULIA_NUM_THREADS=7 julia --project=. -p 4 test_case_1.jl
```

Case 1 on Windows:

```
julia --project=. -p 4 test_case_1.jl
```

Case 3 on macOS:

```
JULIA_NUM_THREADS=7 julia --project=. -p 4 test_case_3.jl
```

Case 3 on Windows:

```
julia --project=. -p 4 test_case_3.jl
```

#### Special test for Case 1

For  Case 1 only, as a test, to ensure that indicator variables and sink flux variables are connected via coupling variables, comment out the following line:

```
optimize_case_1_result = optimize_case_1(case1_ct)
```

uncomment the following lines

```
optimize_case_1_result =
    optimize_case_1(case1_ct; force_first_sink_on = true, force_first_sink_lb = 0.1)
```

and look for the following output (or something similar)

```
Info: Case 1: Optimize constraint tree
objective = 2.000000000000
R_UNKNOWN_SK_DOWN_10fthf_c   flux = 0.100000000000   indicator = 1.000000000000
R_UNKNOWN_SK_UP_10fthf_c   flux = -0.100000000000   indicator = 1.000000000000
FORCED R_UNKNOWN_SK_DOWN_10fthf_c   flux = 0.100000000000   indicator = 1.000000000000
```

Buried somewhere in the middle of the output.
