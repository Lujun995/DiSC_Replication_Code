# DiSC_Replication_Code

This folder contains scripts, data, and R markdown files that can be used to reproduce the results presented in our paper. In the main folder, you will find four markdown files and their corresponding outputs. Additionally, there are several subfolders used to organize data, source code, scripts, parameters, R markdown files, and more.

## Simulations Based on Real-World Data

- **R Markdown:** Real_data_simulation.Rmd
- **R Markdown Output:** Real_data_simulation.html
- **Corresponding Data Folder:** real_data_sim
- **Zipped files:** type1_0715.zip.001 and type1_0715.zip.002 need to be unzipped before use.
- **Description:** These simulations are based on cell type "L2/3" in the prefrontal cortex (PFC) region, sourced from an autism study conducted by Velmeshev, et al. (2019). We permuted labels for the variable `diagnosis` (ASD vs. control) and applied various statistical methods to investigate the association between diagnosis and single-cell RNA sequencing data, with a focus on examining type I error rates and false discovery rates.

## Simulations Based on Parametric Models

- **R Markdown:** Parametric_simulation.Rmd
- **R Markdown Output:** Parametric_simulation.html
- **Corresponding Data Folder:** para_data_sim
- **Zipped files:** results_n12_12_d1_c375_vFALSE_0811.zip.001 and results_n12_12_d1_c375_vFALSE_0811.zip.002 need to be unzipped before use.
- **Description:**
  - The parametric simulation begins with parameters estimated from real-world data, derived from the autism study by Velmeshev et al. (2019) using a deep count autoencoder (DCA) as described by Eraslan et al. (2019) at the cell level. These estimations are then aggregated at the individual level for the simulation.
  - The parametric simulation encompasses four scenarios and three types of changes. The four scenarios include increased effect sizes, increased cell numbers, increased sample sizes, and a no-signal scenario designed to examine type I error and false discovery rate (FDR). The three types of changes are differences in mean, variance, and both mean and variance.
  - The primary objective is to assess the statistical power, type I error rate, and false discovery rate of different methods.

## Real-World Data Analyses

- **R Markdown:** Real_data_analysis.Rmd
- **R Markdown Output:** Real_data_analysis.html
- **Corresponding Data Folder:** real_data_anal
- **Description:** The data used in these analyses come from 17 distinct cell types within the prefrontal cortex region, as part of the autism study by Velmeshev et al. (2019). These analyses aim to test the association between autism spectrum disorder (ASD) and single-cell RNA sequencing data, while adjusting for covariates such as `age`, `sex`, sequence batch (`Seqbatch`), and RNA Integrity Number (`RIN`). The analyses employ various methods, and the results are compared.

## Supplementary Information and Analyses

- **R Markdown:** Supplementary_information.rmd
- **R Markdown Output:** Supplementary_information.html
- **Corresponding Data Folder:** supp_info_anal
- **Description:** The example data set used for these analyses was generated in the "Simulations Based on Parametric Models" step. It includes 8,000 genes and 12 cases and 12 controls, each with 375 cell replicates. The objectives are to (i) determine the number of permutations required for DiSC; (ii) compare the computational time of different methods; (iii) evaluate memory usage by DiSC; and (iv) examine the impact of dropping a proportion of insufficiently sequenced cells during the rarefaction step on the results.

## Other Subdirectories

- **ideas_pipeline-main:** This folder contains scRNA-seq data from the prefrontal cortex (PFC) region of both autism patients and healthy controls, sourced from Velmeshev et al. (2019). The folder also hosts source code necessary for the DCA-IDEAS method, as discussed in Zhang et al. (2022). These contents were adapted from the [ideas_pipeline](https://github.com/Sun-lab/ideas_pipeline/) by Zhang et al. (2022).
- **results:** This directory houses figures and tables generated from the R markdown files above, which are included in our manuscript.
- **source_code:** This folder contains the source code needed to reproduce the tables and figures from the paper. Please note that these files **cannot** be directly executed, as they require additional parameter specifications before running. The complete scripts can be found in the folders mentioned above: real_data_sim/, para_data_sim/, real_data_anal/, and supp_info_anal.

## References

- Eraslan, G., et al. Single-cell RNA-seq denoising using a deep count autoencoder. *Nature communications*, 10.1 (2019): 1-14.
- Velmeshev, D., et al. Single-cell genomics identifies cell type–specific molecular changes in autism. *Science*, 364.6441 (2019): 685-689.
- Zhang, M., et al. IDEAS: individual level differential expression analysis for single-cell RNA-seq data. *Genome Biol*, 23, 33 (2022). [DOI](https://doi.org/10.1186/s13059-022-02605-1).
