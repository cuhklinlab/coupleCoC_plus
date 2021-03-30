# A quick guide to coupleCoC+

## 1. UMAP visualization for the raw data

```R
## UMAP visualization for the datasets in example 2 in the paper
# You need to install the following R packages first
library('R.matlab')
library(umap)
library(RColorBrewer)
source("R/heatmap.R")
colors = brewer.pal(8, "Set3")
# load data
X<-readMat("data/methy_X.mat");X<-as.matrix(X[[1]]);
Y<-readMat("data/methy_Y_link.mat");Y<-as.matrix(Y[[1]]);
Cx_truth<-readMat("data/methy_Cx_truth.mat");Cx_truth<-as.matrix(Cx_truth[[1]]);
Cy_truth<-readMat("data/methy_Cy_truth.mat");Cy_truth<-as.matrix(Cy_truth[[1]]);
# plot source data
Cx_truth[which(Cx_truth==1)]= rep("L4", length(which(Cx_truth==1)))
Cx_truth[which(Cx_truth==2)]= rep("L2/3 IT", length(which(Cx_truth==2)))
colors_CX = colors[1:2]
umap_x = umap(X)
plot.umap(umap_x, labels=as.factor(Cx_truth), main = "A UMAP visualization of scRNA-seq data", colors = colors_CX)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/ex2_S.png "Source data")





```MATLAB
clear
clc
close all

%%load data (take example 3 for example)
load('data/S.mat');load('data/T.mat');load('data/U.mat');
load('data/S_cell_label.mat');load('data/T_cell_label.mat'); 

%%coupleCoC+
%setting hyperparameters
nrowcluster1=2;nrowcluster2=2;ncolcluster=5;ncolcluster0=8;iter=20;
lambda=0.1;beta=0.6;gamma=1;nsub=2;
[Cx, Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0, obj, matm] = coupleCoC_plus(p,q,q0,nrowcluster1,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta,gamma,nsub);

%%results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy);
disp(matm)
```
