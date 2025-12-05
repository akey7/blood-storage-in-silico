# blood-storage-in-silico
Code and data to analyze in silico models of blood storage.

## Prepare to process data

### `input/` folder

The `input/` folder comes preopulated with the following files:

1. [Supplementary data sheet 1 as a `.csv` file from Nemkov et al](https://www.frontiersin.org/api/v4/articles/833242/file/Data_Sheet_1.CSV/833242_supplementary-materials_datasheets_1_csv/1), which is the metabolomics data being analyzed.

2. [RBC-GEM.json v1.3.0 from Haiman et al](https://github.com/z-haiman/RBC-GEM/blob/1.3.0/model/RBC-GEM.json) which is the GEM onto which the metabolomics data above are mapped.

3. `Proportination Sheet 2.csv` which maps columns from the metabolomics data, splits apart columns that contain multiple RBC-GEM metabolites, and proportionates the intensity values among multiple metabolites (if needed), and maps RBC-GEM identifiers to names in the metabolomics data.

4. `Subsystem Category Map.csv`, which maps GEM subsystems into categories for better data visualization. This is the first two columns of [`subsystems.tsv` v1.3.0 of the RBC-GEM](https://github.com/z-haiman/RBC-GEM/blob/1.3.0/data/curation/subsystems.tsv)

### `output/` folder

Please create an `output/` folder.

As configured in the `.gitignore` the contents of this folder will not be synced to GitHub.

Under the `output/` folder, create the following subfolders:

1. `output/c_means_plots/`: This will hold c-means plots.

2. `output/plots/`: This will hold time course plots.

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

## Works Cited

> Haiman, Z. B., Key, A., D’Alessandro, A. & Palsson, B. O. RBC-GEM: A genome-scale metabolic model for systems biology of the human red blood cell. PLoS Comput Biol 21, e1012109 (2025).

> Nemkov, T., Yoshida, T., Nikulina, M. & D’Alessandro, A. High-Throughput Metabolomics Platform for the Rapid Data-Driven Development of Novel Additive Solutions for Blood Storage. Front. Physiol. 13, 833242 (2022).

