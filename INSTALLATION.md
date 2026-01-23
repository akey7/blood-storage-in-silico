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
└── Subsystem Category Map.csv
```
