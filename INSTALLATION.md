# Installation

### Install Julia for Your Platform

[Installation instructions can be found on the language's homepage.](https://julialang.org/)

### Obtain and Install Input Files

1. Obtain `Absolute Quant Data Sheet.xlsx` from the author of this package because it can't be distributed publicly. Place this file in the `input/` folder.

2. Obtain [`RBC-GEM.xml` from the RBC-GEM (Haiman et al)](https://github.com/z-haiman/RBC-GEM/blob/main/model/RBC-GEM.xml) and place in the `input/` folder.

Most of the files are tracked in the git repo. The complete `input/` repo should have this structure.

```
input
├── Absolute Quant Data Sheet.xlsx
├── Data Sheet 1.CSV
├── Proportionation Sheet 2.csv
├── RBC-GEM.json
├── RBC-GEM.xml
├── flux_bounds_overrides.csv
└── Subsystem Category Map.csv
```

Incidentally, here are purposes of select input files:

1. `Data Sheet 1.CSV` which is the metabolomics data being analyzed.

2. `RBC-GEM.*` is the GEM onto which the metabolomics data above are mapped.

3. `Proportination Sheet 2.csv` maps columns from the metabolomics data, splits apart columns that contain multiple RBC-GEM metabolites, and proportionates the intensity values among multiple metabolites (if needed), and maps RBC-GEM identifiers to names in the metabolomics data.

4. `Subsystem Category Map.csv`, maps GEM subsystems into categories for better data visualization. This is the first two columns of [`subsystems.tsv` v1.3.0 of the RBC-GEM](https://github.com/z-haiman/RBC-GEM/blob/1.3.0/data/curation/subsystems.tsv)

5. `flux_bounds_overrides.csv`: Flux bounds in this file override what is specified in the RBC-GEM.

### Create the `output/` Folders

There are a lot of modules and scripts in this repo, and they produce a lot of output files. These files go into the `output/` folder and folders nested within it. Create the `output/` folder and the following subfolders:

```
output
├── pca_plot_dfs
├── pca_plots
├── plots
├── regression_plots
├── relative_absolute_c_means
├── relative_absolute_plots
├── uFBA_histograms_v2
├── gem_dfs
├── ufba_models
```

### Install Dependencies

To install and precompile the Julia dependencies, open a command line in the root of the repo and type the following commands:

```
bash-3.2$ julia --project=.
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.4 (2026-01-06)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |

julia> ]
(BloodStorageInSilico) pkg> instantiate
```

This will instantiate the environment and download the dependencies. After the packages are installed, type backspace and `exit()`.

Further documentation on executing the scripts are found elsewhere in the documentation.

### Build the Documentation

The docstrings are rendered into serachable html pages with a subproject using [Documenter.jl](https://documenter.juliadocs.org/stable/)

```
bash-3.2$ julia --project=docs/
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.4 (2026-01-06)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |

(BloodStorageInSilico/docs) pkg> instantiate
```

After this step is complete build the docs with the following commands from the root of the repo:

```
cd docs/
julia --project=. make.jl
```

When complete, you can open the doucmentation from the following html file relative to the root of the repo:

```
docs/build/index.html
```

Which will present you with nicely formatted docstrings for the functions in the modules.

### Note for Windows

The scripts in this project execute on multiple threads to increase performance. By default, only one thread/core is used. To enable Julia to use all cores in the machine, a reasonable value in the JULIA_NUM_THREADS environment variable must be set. On Windows, you can do this at the user account level in the system settings. For example, on a 64-core machine, you can set JULIA_NUM_THREADS to be the following:

```
JULIA_NUM_THREADS=64
```

On macOS and Linux, you can set the number of threads on the command line, eliminating the need for this extra configuration step. 
