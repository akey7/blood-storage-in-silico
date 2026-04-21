# Blood Storage in Silico - uFBA Analysis Workflow

## UfbaSamplerAnalysisAndViz

Because the `UfbaSampler.jl` module was becoming huge, I split the visualization and analysis functions for `UfbaSampler.jl` into their own module.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSamplerAnalysis]
Order   = [:function]
```

## UfbaSamplerViz3D

3D plots for the uFBA sampler results.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSamplerViz]
Order   = [:function]
```

## ModelGraph

This module contains functions to process metabolic networks as graphs.

```@autodocs
Modules = [BloodStorageInSilico.ModelGraph]
Order   = [:function]
```
