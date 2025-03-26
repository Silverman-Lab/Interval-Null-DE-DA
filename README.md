
## Purpose

This repo contains updated code for BMC Bioinformatics for the paper "Replacing Normalizations with Interval Assumptions Enhances Differential Expression and Differential
Abundance Analyses"

See prior release for code used for the preprint paper.

## Notes on Requirements

This requires the INDExA package version 0.1.0 which can be found at https://github.com/Silverman-Lab/INDExA.
The version is tagged here: https://github.com/Silverman-Lab/INDExA/releases/tag/v0.1.0.

For the CCRCC analysis, the following files must be downloaded and unzipped to enable analysis:

1. [GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE62944&format=file&file=GSE62944%5F01%5F27%5F15%5FTCGA%5F20%5FCancerType%5FSamples%2Etxt%2Egz)
2. [GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE62944&format=file&file=GSE62944%5F06%5F01%5F15%5FTCGA%5F24%5FNormal%5FCancerType%5FSamples%2Etxt%2Egz)
3. [GSM1536837_06_01_15_TCGA_24.tumor_Rsubread_FeatureCounts.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1536837&format=file&file=GSM1536837%5F06%5F01%5F15%5FTCGA%5F24%2Etumor%5FRsubread%5FFeatureCounts%2Etxt%2Egz)
4. [GSM1697009_06_01_15_TCGA_24.normal_Rsubread_FeatureCounts.txt](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM1697009&format=file&file=GSM1697009%5F06%5F01%5F15%5FTCGA%5F24%2Enormal%5FRsubread%5FFeatureCounts%2Etxt%2Egz)

## Notes on Analysis

1. Under `project_code/output` the intermediate files can be found, thus re-analysis is not needed. However, here are details on how to recreate these files and make figures if needed:
   1. figure_2: First run `bash run_all.sh`, then run `Rscript make_plot_2.R`, then make the figure using the ipython notebook
   2. figure_3: First run `Rscript make_plot_3_data.R`, then make the figure using the ipython notebook
   3. figure_s2: First run `Rscript make_plot_S2_data.R`, then run `Rscript make_plot_S2.R`, then make the figure using the ipython notebook

