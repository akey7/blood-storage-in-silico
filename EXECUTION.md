## Purposes of scripts

### `src/` folder

The `src/` folder contains most source code for the analysis that is stored in Julia modules. These modules are not meant to be executed directly, but rather are called by scripts in the root folder. The modules are

1. `src/BloodStorageInSilico.jl`: Enables the `BloodStorageInSilico` module to be compiled during package management.

2. `src/MetaboliteTimelines.jl`: Plots time series of metabolite intensities normalized the the median intensity of that metabolite across all treatments and days.

3. `src/TreatmentsAgainstControlMedians.jl`: Creates timeseries of intensities. Normalization of these intensities is different than in `MetaboliteTimelines.jl`. In `TreatementsAgainstControlMedians.jl`, intensities are normalized by the median of itensity of the "01-Ctrl AS3" intensity for each day. This module performs c-means clustering and visualization of the clustering results.

### Root folder

There are two scripts in the root folder. They are:

1. `plot_metabolite_timelines.jl`: Uses `src/MetaboliteTimelines.jl` to generate plots of timerseries of metabolite timelines.

```
JULIA_NUM_THREADS=7; julia plot_metabolite_timelines.jl
```

2. `treatments_against_control_medians.jl`: Performs c-means clustering with `src/TreatmentsAgainstControlMedians.jl`. On macOS, execute with:

```
JULIA_NUM_THREADS=7; julia treatments_against_control_medians.jl
```

Note that in both macOS commands, the `JULIA_NUM_THREADS` environment variable sets the number of threads that Julia will attempt to use to execute the task. Customize according to your execution environment.

The commands to launch the scripts on Windows are similar, but the JULIA_NUM_THREADS environment variable is not specified in the command line. Rather, the environment variable is configured in settings.

### `ufba_sampler.jl`

To execute on windows (adjsut nchains and concurrent worker processes `-p` according to system architecture). Note that `JULIA_NUM_THREADS` must be set to the appropriate number of threads in settings.

```
julia --project=. -p 7 .\ufba_sampler.jl --nchains 5
```
