# Analysis pipeline files

This repository includes the analysis pipeline of Guo et al., "Single-cell virology.

## `.R` Files

- `generateTidyData.R`: Generate tidy data files from the raw data files
- `sicegarv02_app.R`: This is the pipeline file that runs sicegar package to classify GFP signals as *infection* or *infection and lysis*
- `sicegarv02_app_error.R`: The file for debugging