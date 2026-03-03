# Commands to Execute Scripts

## Multithreaded Execution

Note that in both macOS commands, the `JULIA_NUM_THREADS` environment variable in the execution command sets the number of threads avaiable to Julia.

The commands to launch the scripts on Windows are similar, but the JULIA_NUM_THREADS environment variable is not specified in the command line. Rather, the environment variable is configured in settings as described in [INSTALLATION.md](INSTALLATION.md).

Either way, customize the `JULIA_NUM_THREADS` environment variable to match your environment.

## Setting Up `input/` and `output/` Folders

Before you execute these scripts, the `input/` and `output/` folder must be configured as described in [INSTALLATION.md](INSTALLATION.md).

## uFBA Workflow - Order of Script Execution

The scripts take input and write output files. Some scripts rely on output files previously written by other scripts. The order of script execution presented here maintains the order of reliance of the scripts on each other if such order is important. Such an arrangment of scripts might not be ideal, but it works for now!

## Executing Scripts

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

### (2) `raw_relative_intensities.jl`: 3D PCA Plots of Relative Quant Data

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

On a macOS or Linux machine with 14 cores, an example command to set the number of workers and threads on the same line would be:

```
JULIA_NUM_THREADS=7 julia --project=. -p 7 ufba_sampler.jl --nchains 5 --nmodels -1
```

On a Windows machine with 64 cores, an example to work with your previously set `JULIA_NUM_THREADS` environment variable would be:

```
julia --project=. -p 32 .\ufba_sampler.jl --nchains 5 --nmodels -1
```

Which would sample all models with 10 chains, run all models, and use 32 concurrent workers.

Customize workers, threads, number of chains, and number of models your use case. For quick runs, set the number of models and chains to be small numbers.

Outputs the following files:
1. `outputs/ufba_sampling_status.csv`: That statuses of each uFBA sampling job (fail or ok)
2. `outputs/ufba_sampling.csv`: The samplings of the fluxes. Used by next step.
3. `outputs/case3_sinks.csv`: The unmeasured metabolite sinks from the Case 3 analysis, and whether they had zero or non-zero flux.
4. `outputs/case3_sinks_aggregated.csv`: Counts of zero and nonzero sinks. 
5. `outputs/fba_model_metabolites.csv`: Metabolite ids of the FBA models created for the uFBA runs.

### (5) `ufba_sampler_analysis_and_viz.jl`: Analyze and visualize the results of the uFBA Runs

Outputs csv and plot files of analyses of the uFBA results.

Runs code in the `src/UfbaSamplerAnalysisAndViz.jl`. Shows nifty status bars to indicate progress.

Uses the following input file:

1. Reads the uFBA sampling results file at `output/ufba_sampling.csv`.

Outputs the following files:

1. Writes a reaction id to reaction string yaml file to `output/ufba_sampling.csv`.
2. writes histograms of sampling results (one plot per reaction) to `output/uFBA_histograms_v2/`.
3. Diagnoses the output of the models sampled by uFBA to help find potential problems and writes the diagnostics in `output/ufba_diagnostics.csv`.
4. Writes net fluxes of each pair of sinks to `output/net_sink_fluxes.csv`.
5. Writes a data matrix of median fluxes to `output/flux_vector_data_matrix.csv`.

There are no threads or workers in this script, so execution is simple.

On macOS or Windows:

```
julia --project=. ufba_sampler_analysis_and_viz.jl
```

### (6) `mode_graph.jl`: Analyze the uFBA models as graphs

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

## Other Scripts

There are other scripts that you can run in this project. They are outside of the main uFBA workflow, and are thus optional. They are documented here for completeness.

### `fba_model_builder.jl`: Sanity check on a single FBA model

This script creates a single FBA model and samples the fluxes as a quick sanity check when building models in `FbaModelBuilder.jl`. If it doesn't report any errors, that is good news! However, try the models in the uFBA workflow to check for errors before you run the victory lap.
