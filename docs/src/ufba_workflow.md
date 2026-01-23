# Blood Storage in Silico - uFBA Workflow

uFBA simulation is the main workflow for this system. The modules in this section run this workflow. Other modules were written for more exploratory purposes outside of this workflow.

I present the modules in the order they are used in the workflow.

## RawRelativeIntensities

RawRelativeIntensities works with the relative quantification data. This includes exploratory analysis of the relative intensity data and preparation to incorporate with absolute quantification data.

The paper it is designed to work with is from the following paper:

> Nemkov, T., Yoshida, T., Nikulina, M. & D’Alessandro, A. High-Throughput Metabolomics Platform for the Rapid Data-Driven Development of Novel Additive Solutions for Blood Storage. Front. Physiol. 13, 833242 (2022).

```@autodocs
Modules = [BloodStorageInSilico.RawRelativeIntensities]
Order   = [:function]
```

## AbsoluteQuant

The original blood storage study by Nemkov et al (2022) used relative quantification for the study. However, uFBA analyses require absolute quantification. In the absence of absolute quant data, we had to make an approximation of absolute quant values. This approximation is made by combining the relative quant data with absolute quant data of blood initially stored in similar conditions (1 week, AS3 additive solution, just like the control of the realtive quant study). Using the aboslute quant numbers as a baseline, we scaled the fold changes relative quantification study byt the absolute values to approximate an absolute quant study of blood storage metabolites over the time course of the relative quantification study.

Also, in mass spec, sometimes multiple compounds appear as a single peak and cannot be resolved. However, as required by the uFBA modeling, these peaks must be seprated into single compounds. This is handled by a "proporination" file that serves two functions: (1) it maps compound names to metabolite_ids as used in the RBC-GEM by Haiman et al and (2) seprates combined peaks into individual compounds, each with a proportional fraction of the absolute concentration.

With all these assumptions taken together, `AbsoluteQuant.jl` performs the calculations to make these approximations to feed into the uFBA study

```@autodocs
Modules = [BloodStorageInSilico.AbsoluteQuant]
Order   = [:function]
```

## UfbaSampler

UfbaSampler creates a three-pathway model and samples it for an unsteady flux balance analysis (uFBA) study as described by Brodbar et al.

> Bordbar, A. et al. Elucidating dynamic metabolic physiology through network integration of quantitative time-course metabolomics. Sci Rep 7, 46249 (2017).

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler]
Order   = [:function]
```

## UfbaSamplerViz

Because the `UfbaSampler.jl` module was becoming huge, I split the visualization functions for `UfbaSampler.jl` into their own module.

```@autodocs
Modules = [BloodStorageInSilico.UfbaSamplerViz]
Order   = [:function]
```