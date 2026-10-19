# XRF and PCA Analysis of Siliceous Raw Materials (La Sima Cave)

This repository contains the R code and geochemical data (XRF) used to assess the compositional affinity between Middle Palaeolithic lithic artefacts from La Sima Cave (Constantina, Seville, Spain) and the geological quartz sources in its surrounding landscape.

## Description

The script performs a multivariate statistical analysis on X-ray fluorescence (XRF) data expressed as weight percentages of oxides. Given the compositional nature of the data, a centred log-ratio (CLR) transformation and Z-scoring are applied prior to Principal Component Analysis (PCA). 

Two approaches are executed:
1. A PCA including all elements (including P₂O₅ and CaO).
2. A PCA excluding P₂O₅ and CaO, to control for the effect of post-depositional contamination (anthropised Terra Rossa).

Finally, each artefact is assigned to its nearest geological source using Euclidean distance in the space of the first five principal components.

## Requirements

To run the code, **R** (version 4.0 or higher) must be installed. The script automatically installs and loads the following packages:
* `ggplot2`
* `dplyr`
* `tidyr`
* `compositions`
* `ggrepel`
* `RColorBrewer`
* `pheatmap`
* `psych`

## Usage

1. Download the `Script_XRF_Analysis.R` file.
2. Open RStudio or your preferred R environment.
3. Run the entire script. 
4. The script will automatically generate an output folder on your desktop (or working directory) containing all high-resolution figures (300 dpi) and CSV files with the assignment results.

## Data

The input data is embedded directly in the code to ensure reproducibility. It comprises 13 samples: 4 geological sources (CZ) and 9 archaeological artefacts (SIM). CO₂ values (balance of light elements) were excluded from all statistical analyses.

## Citation

If you use this code or the derived results, please cite:

> Medina Luque, E. J. (2026). *XRF and PCA Analysis of Siliceous Raw Materials (La Sima Cave)* (v1.0.0). Zenodo.

## Author

**Eusebio Jesús Medina Luque**
University of Córdoba
ORCID: [0000-0003-4986-7089](https://orcid.org/0000-0003-4986-7089)