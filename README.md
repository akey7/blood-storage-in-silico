# blood-storage-in-silico
Code and data to analyze in silico models of blood storage.

## Prepare to process data

## `input/` folder

The `input/` folder comes preopulated with the following files:

1. [Supplementary data sheet 1 as a `.csv` file from Nemkov et al](https://www.frontiersin.org/api/v4/articles/833242/file/Data_Sheet_1.CSV/833242_supplementary-materials_datasheets_1_csv/1), which is the metabolomics data being analyzed.

2. [RBC-GEM.json v1.3.0 from Haiman et al](https://github.com/z-haiman/RBC-GEM/blob/1.3.0/model/RBC-GEM.json) which is the GEM onto which the metabolomics data above are mapped.

3. `Proportination Sheet 2.csv` which maps columns from the metabolomics data, splits apart columns that contain multiple RBC-GEM metabolites, and proportionates the intensity values among multiple metabolites (if needed), and maps RBC-GEM identifiers to names in the metabolomics data.

4. `Subsystem Category Map.csv`, which maps GEM subsystems into categories for better data visualization. This is the first two columns of [`subsystems.tsv` v1.3.0 of the RBC-GEM](https://github.com/z-haiman/RBC-GEM/blob/1.3.0/data/curation/subsystems.tsv)

## Works Cited

> Haiman, Z. B., Key, A., D’Alessandro, A. & Palsson, B. O. RBC-GEM: A genome-scale metabolic model for systems biology of the human red blood cell. PLoS Comput Biol 21, e1012109 (2025).

> Nemkov, T., Yoshida, T., Nikulina, M. & D’Alessandro, A. High-Throughput Metabolomics Platform for the Rapid Data-Driven Development of Novel Additive Solutions for Blood Storage. Front. Physiol. 13, 833242 (2022).

