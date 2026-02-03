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

Uses `src/RawRelativeIntensities.jl` to make 3D PCA plots reducing relative metabolite abundances down to fewer features.

This script does not use multiple workers or threads, so executing it is easy.

On macOS or Windows:

```
julia --project=. raw_relative_intensities.jl
```

This will display interactive GLMakie scatter plots of the first 3 principal components. Screen capture to obtain files for publication or presentations.

### (3) `absolute_quant.jl`: Approximate Absolute Quantifications and Regress Concentration Change Rates

Uses `src/AbsoluteQuant.jl` to perform the following tasks:

1. Loads data relative and absolute quant data files.

2. Writes quality checks to `output/qc_fold_changes.csv` and `output/qc_fold_change_zeros.csv`.

3. Map and proportinoate (see documentation built in [INSTALLATION.md](INSTALLATION.md)) combined metaboilite names to single metabolite ids from the RBC-GEM. Approximate absolute concentrations for the relative quant data using accompanying absolute quant data. Write result to `output/relative_absolute_quant.csv`.

4. Plots the approximated absolute concentrations over time to plots in the `output/relative_absolute_plots/[cleaned metabolite name].png`

5. Perform c-means clustering on the metabolite trajectories. Create plots of c-means clusters and an accompanying elbow plot for each additive to `output/relative_absolute_c_means/`. Write cluster memberships to `output/c_means_primary_clusters.csv`.

6. Performs PCA analysis on the approximate absolute quant values **NOTE: This functionality is deprecated in preference of the PCA in `raw_relative_intensities.jl` file**

7. Performs the regression to determine the rates of metabolite concentration changes and writes the result to `output/concentration_rates.csv`. Plots regression data and stores the plots in `output/regression_plots/`

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
JULIA_NUM_THREADS=7 julia --project=. -p 7 .\ufba_sampler.jl --nchains 10 --nmodels -1
```

On a Windows machine with 64 cores, an example to work with your previously set `JULIA_NUM_THREADS` environment variable would be:

```
julia --project=. -p 32 .\ufba_sampler.jl --nchains 10 --nmodels -1
```

Which would sample all models with 10 chains, run all models, and use 32 concurrent workers.

Customize workers, threads, number of chains, and number of models your use case. For quick runs, set the number of models and chains to be small numbers.

### (5) `ufba_sampler_analysis_and_viz.jl`: Analyze and visualize the results of the uFBA Runs

Runs code in the `src/UfbaSamplerAnalysisAndViz.jl`. Reads the uFBA sampling results file at `output/ufba_sampling.csv`, writes a reaction id to reaction string yaml file to `output/ufba_sampling.csv`, and writes histograms of sampling results (one plot per reaction) to `output/uFBA_histograms_v2/`. Makes a nifty progress bar to show progress. Also diagnoses the output of the models sampled by uFBA to find potential problems. 

There are no fancy threads or workers here, so execution is simple.

On macOS or Windows:

```
julia --project=. ufba_sampler_viz.jl
```
