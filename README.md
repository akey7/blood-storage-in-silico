# Blood Storage *In Silico*

This project studies refrigerated blood storage with unsteady flux balance analysis (uFBA).

## Documentation

1. [Installation instructions](INSTALLATION.md), covering installing dependencies and formatting docstrings.

2. [Explanations of an execution instructions](EXECUTION.md) for the scripts.

3. [Works cited that are important to this project](WORKS_CITED.md)

4. [New to Julia documentation](NEW_TO_JULIA.md) for users new to the Julia programming language.

5. [Guide for reviewer](GUIDE_FOR_REVIEWERS.md) for a guide on checking the code for the purposes of scientific peer review.

6. **Once built**, `docs/build/index.html` contains documentation about the modules in `src/` and functions in those modules that encode the functionality of this repo.

## Folder Structure

Relative to the root of the repo, the top level folders and groups of files are listed below along with breif explanations of their purpose.

1. `docs/`: Subproject to build and hold nicely-formatted versions for the docstrings for the functions in the modules. Built versions of the documentation are not contained in this repo by default. See [INSTALLATION.md](INSTALLATION.md) for instructions on building the documentation.

2. `input/`: Holds input data files for that the scripts read during operation. Not all input files can be committed to this repo due to intellectual property concerns. For more information in completing this part of the repo, see [INSTALLATION.md](INSTALLATION.md).

3. `output/`: A directory tree where script output is stored. There are number of subfolders that need to exist for all the data to be written. See setup information in [INSTALLATION.md](INSTALLATION.md) for instructions on setting up this folder.

4. `src/`: Modules that enable the scripts to run. These modules are not meant to be executed directly; rather, they are loaded as dependencies of other scripts that are meant to be executed from the command line. For more information on the modules themselves, build the documentation as instructed in [INSTALLATION.md](INSTALLATION.md).

5. `*.jl`: Command line tools meant to execute the modules in this project to read input files, write output files, and generate research results. See more information about how to execute these files properly in [EXECUTION.md](EXECUTION.md).

6. `*.md`: Documentation files.

7. `Project.toml`: Julia project file that tracks dependencies.
