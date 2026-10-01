# Blood Storage in Silico - Quantification Workflow

The quantification workflow does an exploratory data analysis (EDA) of the relative and absolute qauntification metabolomics data. It also combines these data sources to approximate absolute quantification which is needed by the uFBA workflow.

They are presented here in the order they are executed to produce results.

## MetaboliteTimelines

This module contains functionality to visualize the median of relative quant data intensities and find correlations between metabolites. The syntax of the DataFrame manipulations in this module is not as eloquent as it is in other modules, but it gets the job done.

```@autodocs
Modules = [BloodStorageInSilico.MetaboliteTimelines]
Order   = [:function]
```

## RawRelativeIntensities

RawRelativeIntensities works with the relative quantification data. This includes exploratory analysis of the relative intensity data and preparation to incorporate with absolute quantification data.

The relative quantification metabolomics it is designed to work with is from the following paper:

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
