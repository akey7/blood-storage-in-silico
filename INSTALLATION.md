# Installation

## Install Julia for Your Platform

[Installation instructions can be found on the language's homepage.](https://julialang.org/downloads/)

## Configure `input/` folder with necessary data

There are many data required in `input/` to run this code. The maximum number of publicly accessible input data files are committed to this repo directly. However, **there are proprietary input files that cannot be committed to a public repo and must be obtained from the paper authors separately.** Files that must be obtained from the authors are listed as "distributed to trusted reviewers/collaborators" in the table below.

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

## Install Dependencies

To install and precompile the Julia dependencies, open a command line in the root of the repo and type the following commands **(Note that you need to type `]` at the `julia>` prompt as shown below)**:

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

This will instantiate the environment, download the dependencies, and compile them. **Compilation takes approximately 15 minutes.** After the packages are installed, type backspace and `exit()`.

Further documentation on executing the scripts are found elsewhere in [EXECUTION.md](EXECUTION.md).

## Build the Documentation

The docstrings are rendered into serachable html pages with a subproject using [Documenter.jl](https://documenter.juliadocs.org/stable/). These documentation pages need to be built during installation.

First, launch a Julia REPL in the `docs/` project as shown below.

```
bash-3.2$ julia --project=docs/
```

 **Note that you need to type `]` at the `julia>` prompt as shown below**

```
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.12.4 (2026-01-06)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org release
|__/                   |

julia> ]
(BloodStorageInSilico/docs) pkg> instantiate
```

This will instantiate the environment and download the dependencies **to build the documentation**. After the packages are installed, type backspace and `exit()`.

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

## *Windows Only*: Setting a required user environment variable on Windows

The scripts in this project execute on multiple threads to increase performance. By default, only one thread/core is used. To enable Julia to use all cores in the machine, a reasonable value in the JULIA_NUM_THREADS environment variable must be set. On Windows, you can do this at the user account level in the system settings. 

1. Press `Windows + I` to open **Settings**.  
2. In the left sidebar, click **System**.  
3. On the right, scroll down and click **About**.  
4. In the About page, under **Related links**, click **Advanced system settings**. This opens the **System Properties** window.
5. In the System Properties window, make sure the **Advanced** tab is selected.
6. Click the **Environment Variables…** button near the bottom. This opens the **Environment Variables** dialog.
7. In the top section labeled **User variables for \<your username\>**, click **New…**. This ensures the variable is created only for your user account (not system-wide).
8. In **Variable name**, type `JULIA_NUM_THREADS`  
9. In **Variable value**, type the number of threads you would like to allocate. For example on my 64-core AMD Threadripper, I set this value to **64**.
10. Click **OK** to close the New Variable dialog.  
11. Click **OK** to close the Environment Variables dialog.  
12. Click **OK** (or **Apply**) to close System Properties.  
13. Close Settings, then restart any open apps or terminals that need to read the new variable.

You can introduce this in your docs as “Follow these steps on Windows 11 to create a user environment variable” and then list the steps. How would you like to phrase step 7 to make it extra clear that they must use the **User variables** section and not the **System variables** section?

**Ensure that you close and reopen any PowerShell or command prompt windows to ensure these settings take effect.**

On macOS and Linux, you can set the number of threads on the command line, eliminating the need for this extra configuration step. 
