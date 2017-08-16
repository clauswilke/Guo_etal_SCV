# Data and code for Guo et al., Single-cell virology

This repository includes the datsets and analysis pipeline of Guo et al., "Single-cell virology: on-chip investigation of viral infection dynamics". 

The file `SCV_data_information_short.xlsx` includes detailed information about the samples including the conditions of the experiment and the threshold intensity value used in the experiments.

The file `data_figure_relations.txt` includes information that links datasets with figures.

## Folders

- `codeSicegar`: analysis pipeline files, built on top of the R package [sicegar](https://CRAN.R-project.org/package=sicegar)
- `rawData`: raw intensity measurements from single-cell infections
- `tidyData`: same raw intensity files in tidy form
- `detailedSicegarResults`: all the outputs of the sicegar package for each sample set
- `finalResultsRMarkDown`: various files: html files that show comparisons of multiple datasets, files that includes all the outputs of sicegar package combined with experimental condition details, and summary tables that summarize sicegar outputs 
