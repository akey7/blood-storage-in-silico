# Blood Storage in Silico - uFBA Analysis Workflow

## UfbaSamplerAnalysis

Because the `UfbaSampler.jl` module was becoming huge, I split the visualization and analysis functions for `UfbaSampler.jl` into their own module.

Running this module is very time-consuming, so it was split from the second analysis module.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSamplerAnalysis]
Order   = [:function]
```

## UfbaSamplerAnalysis2

3D plots for the uFBA sampler results.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSamplerAnalysis2]
Order   = [:function]
```

## ModelGraph

This module contains functions to process metabolic networks as graphs.

```@autodocs
Modules = [BloodStorageInSilico.ModelGraph]
Order   = [:function]
```
