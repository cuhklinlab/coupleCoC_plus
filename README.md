# The implementation of coupleCoC+
This package contains the source code for the production of results in the paper "coupleCoC+: a information-theoretic co-clustering-based transfer learning framework for the integrative clustering of single-cell genomic data" by Pengcheng Zeng and Zhixiang Lin. The code include matlab code for the coupleCoC+ algorithm and R code for the downstream analysis. 

## Requirements
* MATLAB R2019b for the coupleCoC+ model
- R version 4.0.4 for data visualization and downstream analysis 

## Main functions
Main functions in this package include:
* coupleCoC_plus.m: implementation of the coupleCoC+ algorithm
* swap_label_plus.m: determination of the matched labels between the source data and target data by our criterion
* clu_eval.m: calculate and show the ARI and NMI values of clustering results
* heatmap.R: display the heatmaps of data after clustering by coupleCoC+

## Installation of packages for benchmark methods
To reproduce the results by benchmark methods, you should also install R packages as follows:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SC3", "SIMLR", "BPRMeth")) # packages for methods SC3, SIMLR and BPRMeth
```
```R
install.packages("devtools")#only if you have not installed the package "devtools"
library(devtools)
install_github("shibiaowan/SHARP") # package for method SHARP
install_github("yycunc/SAMEclustering") # package for method SAMEclustering
install_github("aertslab/RcisTarget")
install_github("aertslab/AUCell")
install_github("aertslab/cisTopic") # packages for method cisTopic
install_github("cuhklinlab/scACE") # package for method scACE
```
```R
install.packages(c("rliger","Seurat","data.table","Matrix","proxy","Rtsne","densityClust","data.table","irlba","umap","ggplot2")) # packages for methods LIGER, Seurat, Cusanovich2018 and data visualization
```
The matlab code for the remaining benchmark methods (coupleCoC and CoC) are included in this package.

## Datasets and examples
Please check the [vigenette](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/vignettes/vignette.md) for a tutorial. A few examples are contained for a quick start of coupleCoC+. The analytical scripts for all benchmark methods are also included.

## Acknowledgement
We utilized and modified part of the code of STC algorithm, therefore please refer to the STC code at https://github.com/sykim122/STC.

