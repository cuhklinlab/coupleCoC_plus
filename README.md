# The implementation of coupleCoC+
The source code for the production of results in the paper "coupleCoC+: a information-theoretic co-clustering-based transfer learning framework for the integrative clustering of single-cell genomic data" by Pengcheng Zeng and Zhixiang Lin. The code include matlab code for the coupleCoC+ algorithm and R code for the downstream analysis. 

## Requirements
* MATLAB R2019b for the coupleCoC+ model (including the matlab-based benchmark methods coupleCoC and CoC)
- R version 4.0.4 for data visualization and downstream analysis (including the R-based benchmark methods LIGER, Seurat, SC3, SIMLR, etc.)

## Installation
To reproduce the results by benchmark methods, you should also install R packages as follows:
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c('SC3', "SIMLR")) 
```
```R
install.packages("devtools")#only if you have not installed the package "devtools"
library(devtools)
install_github("shibiaowan/SHARP")
install_github("yycunc/SAMEclustering")
install_github("aertslab/RcisTarget")
install_github("aertslab/AUCell")
install_github("aertslab/cisTopic")
```
```R
install.packages(c('rliger','Seurat'))
```
## Acknowledgement
Because it draws on the code of STC algorithm, please also refer to the STC code at https://github.com/sykim122/STC. Note that because the coupleCoC_plus algorithm requires the random initialization, the clustering results might not be exactly the same as those presented in the paper.
(1) Please run "real_data_results.m" to show the results by coupleCoC+ algorithm for the third real example in the paper.
