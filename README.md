# Basketwork Eel Population Genomics Analysis

Population genetics analysis of Diastobranchus using ddRAD-seq SNP data.

## Contents

- `Analysis_R_code_Basketwork_Eel_ddRAD.R` — main R analysis script for SNP filtering, PCA, AMOVA, FST, and IBD analyses

## Overview

This repository contains the analysis workflow used for the Diastobranchus ddRAD population genomics project. The script performs:

- SNP data import and genlight conversion
- Metadata matching and population assignment
- Filtering for locus call rate and HWE
- PCA and clustering analyses
- AMOVA and pairwise FST calculations
- Isolation by distance analyses
- Mapping and plotting of sample locations

## Requirements

This analysis requires the following R packages:

- adegenet
- vcfR
- SNPRelate
- HardyWeinberg
- dartR.base
- dartR.popgen
- StAMPP
- vegan
- geosphere
- hierfstat
- rnaturalearth
- sf
- ggplot2
- maps

## Input files

The script expects the following files to be placed in the working directory:

- `NXGSQCAGRF25040091-3_variants.vcf`
- `EelGenoInfo2026.csv`
These files are available in the CSIRO data archive: https://data.csiro.au/collection/csiro:78574 

## Usage

1. Update the `Folder` path in the script to your local working directory.
2. Open the script in RStudio or another R environment.
3. Run the script.

## Raw file link for static repositories

Use this link when you want a “latest version” link:

https://raw.githubusercontent.com/bdeagle-code/Basketwork-eel-population-genomics-analysis/main/Analysis_R_code_Basketwork_Eel_ddRAD.R

If you want a fixed version, create a GitHub release or tag and use a versioned URL such as:

https://raw.githubusercontent.com/bdeagle-code/Basketwork-eel-population-genomics-analysis/v1.0/Analysis_R_code_Basketwork_Eel_ddRAD.R

## Notes

This script was used for the initial analysis for a paper submission. Check for any updates/revised code.

## License

This repository is intended for research and reproducibility purposes.
