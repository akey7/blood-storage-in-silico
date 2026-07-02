# Installation

### Install Julia for Your Platform

[Installation instructions can be found on the language's homepage.](https://julialang.org/)

### Obtain and Install Input Files

| File                                          | Required for full reproduction? |                   Public? | How obtained                                   | Used by                          |
| --------------------------------------------- | ------------------------------: | ------------------------: | ---------------------------------------------- | -------------------------------- |
| `Absolute Quant Data Sheet.xlsx`              |                             Yes |                        No | distributed to trusted reviewers/collaborators | `absolute_quant.jl`              |
| `Absolute Quant Extracellular Datasheet.xlsx` |                             Yes |                        No | distributed to trusted reviewers/collaborators | `absolute_quant.jl`              |
| `Data Sheet 1.CSV`                            |                             Yes |                   Yes/repo-tracked | Included in repo                         | `absolute_quant.jl` |
| `AS Dev Library Trial 1.csv`                  |                             Yes |                   No | distributed to trusted reviewers/collaborators                             | `absolute_quant.jl`        |
| `RBC-GEM.xml`                                 |                             Yes |                       Yes/repo-tracked | Included in repo                             | model construction               |
| `flux_bounds_overrides.csv`                   |                             Yes |          Yes/repo-tracked | included in repo                               | FBA model construction           |
| `metabolite_measurement_opt_outs.csv`         |             Yes | Yes/repo-tracked or local | included/template                              | uFBA model construction          |
| `sink_opt_ins.csv`                            |             Yes | Yes/repo-tracked or local | included/template                              | pruning override                 |

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
├── uFBA_densities
├── uFBA_heatmaps
├── gem_dfs
├── ufba_models
├── flux_vector_data_matrices
├── relative_quant_2
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
