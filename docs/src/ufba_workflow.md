# Blood Storage in Silico - uFBA Workflow

uFBA simulation is the main workflow for this system. The modules in this section run this workflow. Other modules were written for more exploratory purposes outside of this workflow.

I present the modules in the order they are used in the workflow.

## FbaModelBuilder

FbaModelBuilder has functions for loading the RBC-GEM and creates the models for uFBA analysis.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler.FbaModelBuilder]
Order = [:function]
```

Many of the reactions in this model are from the following paper:

> Bordbar, A. et al. Identified metabolic signature for assessing red blood cell unit quality is associated with endothelial damage markers and clinical outcomes. Transfusion 56, 852–862 (2016).

Reaction ids have been mapped from that paper, released in 2016, to reaction ids in the RBC-GEM released in 2024.

## MetaboliteBounds

MetaboliteBounds is for configuring metaolite bounds in models before optimizations.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler.MetaboliteBounds]
Order   = [:function]
```

## UfbaSampler

UfbaSampler creates a three-pathway model and samples it for an unsteady flux balance analysis (uFBA) study as described by Brodbar et al.

> Bordbar, A. et al. Elucidating dynamic metabolic physiology through network integration of quantitative time-course metabolomics. Sci Rep 7, 46249 (2017).

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler]
Order   = [:function]
```

## PruningOptimizations

PruningOptimizations performs optimizations on uFBA models to prune unneeded sinks.

> Bordbar, A. et al. Elucidating dynamic metabolic physiology through network integration of quantitative time-course metabolomics. Sci Rep 7, 46249 (2017).

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler.PruningOptimizations]
Order   = [:function]
```
