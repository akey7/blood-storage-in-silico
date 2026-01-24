# Commands to Execute Scripts

## Multithreaded Execution

Note that in both macOS commands, the `JULIA_NUM_THREADS` environment variable in the execution command sets the number of threads avaiable to Julia.

The commands to launch the scripts on Windows are similar, but the JULIA_NUM_THREADS environment variable is not specified in the command line. Rather, the environment variable is configured in settings as described in [INSTALLATION.md](INSTALLATION.md).

Either way, customize the `JULIA_NUM_THREADS` environment variable to match your environment.

## Setting Up `input/` and `output/` Folders

Before you execute these scripts, the `input/` and `output/` folder must be configured as described in [INSTALLATION.md](INSTALLATION.md).

## Executing Scripts

All commands are issued from the root of the repo.

### `plot_metabolite_timelines.jl`: Plot Timelines of Relative Metabolite Intensities

1. `plot_metabolite_timelines.jl`: Uses `src/MetaboliteTimelines.jl` to generate plots of timerseries of metabolite timelines.

on macOS

```
JULIA_NUM_THREADS=7 julia --project=. plot_metabolite_timelines.jl
```

on Windows (assuming your `JULIA_NUM_THREADS` environment variable is set)

```
julia --project=. plot_metabolite_timelines.jl
```

Output will be saved to `output/normalized_abundance_correlations.csv` and `output/plots`

### `ufba_sampler.jl`: Run uFBA Sampling Jobs

In addition to multithrading, the uFBA sampling module uses concurrent worker processes to fully utilize the hardware executing the script. There is an optimal point to set the number of workers: if there are too few workers, the job will take a needlessly long time to execute. With too many workers, the script takes too long to launch.

The command line arguments to the Julia environment and script are the following:

1. `-p`: How many workers Julia will lauch for the sampling task.

2. `--nchains`: The number of chains of samples for each uFBA model. More chains means more sampling and longer execution.

3. `--nmodels`: Number of uFBA models to analyze (-1 for all possible models).

On a macOS or Linux machine with 14 cores, an example command to set the number of workers and threads on the same line would be:

```
JULIA_NUM_THREADS=7 julia --project=. -p 7 .\ufba_sampler.jl --nchains 10 --nmodels 5
```

On a Windows machine with 64 cores, an example to work with your previously set `JULIA_NUM_THREADS` environment variable would be:

```
julia --project=. -p 32 .\ufba_sampler.jl --nchains 10 --nmodels -1
```

Which would sample all models with 10 chains, run all models, and use 32 concurrent workers.

Customize workers, threads, number of chains, and number of models your use case. For quick runs, set the number of models and chains to be small numbers.
