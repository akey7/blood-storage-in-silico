# Installation

### Install Julia for Your Platform

[Installation instructions can be found on the language's homepage.](https://julialang.org/)

### Configure `input/` folder with necessary data

The maximum number of publicly accessible input data files are committed to this repo directly. However, **there are proprietary input files that cannot be committed to a public repo and must be obtained from the paper authors seprately.** Files that must be obtained from the authors are listed as "distributed to trusted reviewers/collaborators" in the table below.

| How obtained                                   | File                                          | Required for full reproduction? |                   Public? | Used by                          |
| ----------------------------------------------- | --------------------------------------------- | ------------------------------: | ------------------------: | -------------------------------- |
| distributed to trusted reviewers/collaborators | `Absolute Quant Data Sheet.xlsx`              |                             Yes |                        No | `absolute_quant.jl`              |
| distributed to trusted reviewers/collaborators | `Absolute Quant Extracellular Datasheet.xlsx` |                             Yes |                        No | `absolute_quant.jl`              |
| Included in repo                         | `Data Sheet 1.CSV`                            |                             Yes |                   Yes/repo-tracked | `absolute_quant.jl` |
| distributed to trusted reviewers/collaborators                             | `AS Dev Library Trial 1.csv`                  |                             Yes |                   No | `absolute_quant.jl`        |
| Included in repo                             | `RBC-GEM.xml`                                 |                             Yes |                       Yes/repo-tracked | model construction               |
| included in repo                               | `flux_bounds_overrides.csv`                   |                             Yes |          Yes/repo-tracked | FBA model construction           |
| included/template                              | `metabolite_measurement_opt_outs.csv`         |             Yes | Yes/repo-tracked or local | uFBA model construction          |
| included/template                              | `sink_opt_ins.csv`                            |             Yes | Yes/repo-tracked or local | pruning override                 |

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
