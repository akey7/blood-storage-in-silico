# Blood Storage in Silico

## Overview

### Additives and Final Time Points

Additives are repeatedly referenced throughout this documentation. These are the additives from the study by Nemkov et al:

> Nemkov, T., Yoshida, T., Nikulina, M. & D’Alessandro, A. High-Throughput Metabolomics Platform for the Rapid Data-Driven Development of Novel Additive Solutions for Blood Storage. Front. Physiol. 13, 833242 (2022).

They are:
1. 01-Ctrl AS3
2. 02-Adenosine
3. 03-Glutamine
4. 04-Methionine
5. 07-NAC
6. 08-Taurine

Metabolites were measured each week of blood storage. In the `AbsoluteQuant.jl` module below, the rates of concentration change were calculated between adjacent weeks, which gives the following `final_time` points as:

- 2: Weeks 1 to 2
- 3: Weeks 2 to 3
- 4: Weeks 3 to 4
- 5: Weeks 4 to 5
- 6: Weeks 5 to 6

## RawRelativeIntensities

RawRelativeIntensities works with the relative quantification data. This includes exploratory analysis of the relative intensity data and preparation to incorporate with absolute quantification data.

The paper it is designed to work with is from the following paper:

> Nemkov, T., Yoshida, T., Nikulina, M. & D’Alessandro, A. High-Throughput Metabolomics Platform for the Rapid Data-Driven Development of Novel Additive Solutions for Blood Storage. Front. Physiol. 13, 833242 (2022).

```@autodocs
Modules = [BloodStorageInSilico.RawRelativeIntensities]
Order   = [:function]
```

## UfbaSampler

UfbaSampler creates a three-pathway model and samples it for an unsteady flux balance analysis (uFBA) study as described by Brodbar et al.

> Bordbar, A. et al. Elucidating dynamic metabolic physiology through network integration of quantitative time-course metabolomics. Sci Rep 7, 46249 (2017).

```@autodocs
Modules = [BloodStorageInSilico.UfbaSampler]
Order   = [:function]
```
