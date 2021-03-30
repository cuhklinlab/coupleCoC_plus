# A quick guide to coupleCoC+

## 1. UMAP visualization for the raw data

```R
## UMAP visualization for the datasets in example 3 in the paper
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
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/ex3_S.png "Source data")

```R
# target data
Cy_truth[which(Cy_truth==1)]= rep("L4", length(which(Cy_truth==1)))
Cy_truth[which(Cy_truth==2)]= rep("L2/3 IT", length(which(Cy_truth==2)))
colors_CY = colors[1:2]
umap_y = umap(Y)
plot.umap(umap_y, labels=as.factor(Cy_truth), main = "A UMAP visualization of sc-methylation data", colors = colors_CY)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/ex3_T.png "Source data")

## 2. Implementation of coupleCoC+

```MATLAB
clear
clc
close all

%%load the processed data in example 3 in the paper
load('data/S.mat');load('data/T.mat');load('data/U.mat');
load('data/S_cell_label.mat');load('data/T_cell_label.mat'); 

%%coupleCoC+
%setting the values of hyperparameters
nrowcluster1=2;nrowcluster2=2;ncolcluster=5;ncolcluster0=8;iter=20;
lambda=0.1;beta=0.6;gamma=1;nsub=2;
[Cx, Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0, obj, matm] = coupleCoC_plus(p,q,q0,nrowcluster1,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta,gamma,nsub);

%% results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy); %% Note that this function produces the contingency table and four metrics values, two of which, i.e. ARI and NMI are utilized in this paper
disp(matm) %% show the matched order of labels in two datasets for each iteration. For example, [1 2 2 1] means the first and the second kind of cell types in the source dataset are matched with the the second and the first kind of cell types in the target dataset, respectively.
```

```MATLAB
>> [TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy) 
TAB_X =
    0    1401
   974    8

TAB_Y =
   386    26
    11   679

Eval_tab =
  4×2 table
                 X          Y   
              _______    _______

    Purity    0.99664    0.96642
    RI        0.99331    0.93505
    ARI        0.9866    0.86946
    NMI       0.97031    0.78253

>> disp(matm) 
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
     1     2     2     1
```
```MATLAB
%%determination of N_sub
gvalue = [];
for N_sub=1:2
[~, ~, dis] = swap_label_plus(cluster_p,cluster_q,N_sub);
gvalue(N_sub)= dis/(N_sub*log(N_sub+1));
end
gvalue %% We set the number of matched N_sub as 2, because the values of the objective function g(N_sub) for choosing N_sub are smaller when N_sub = 2 (0.138 when N_sub = 1 and 0.061 when N_sub = 2.
respectively)
```
```MATLAB
>> gvalue

gvalue =

    0.1382    0.0613
```

## 3. Heatmaps of clustering results by coupleCoC+

```R
## import libraries
source("R/heatmap.R")
library('gplots')
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
library('R.matlab')

## load the raw data and clustering results by coupleCoC+
X0<-readMat("data/methy_X.mat");X0<-as.matrix(X0[[1]]);
Y_link<-readMat("data/methy_Y_link.mat");Y_link<-as.matrix(Y_link[[1]]);
Y_unlink<-readMat("data/methy_Y_unlink.mat");Y_unlink<-as.matrix(Y_unlink[[1]]);

CX_best<-readMat("data/methy_Cx.mat");CX_best<-as.vector(CX_best[[1]])
CY_best<-readMat("data/methy_Cy.mat");CY_best<-as.vector(CY_best[[1]])
CZ_best<-readMat("data/methy_Cz.mat");CZ_best<-as.vector(CZ_best[[1]])
CZ0_best<-readMat("data/methy_Cz0.mat");CZ0_best<-as.vector(CZ0_best[[1]])

## switch label based on the result of matm
old_label_x=c(1,2)
CX_best<-reorder_label(CX_best,old_label_x)
old_label_y=c(2,1)
CY_best<-reorder_label(CY_best,old_label_y)
X_raw = X0;Y_raw_link = Y_link;Y_raw_unlink = Y_unlink;

## obtain the separation lines in the heatmap
gene_order_Z <- clu_sep(CZ_best)
gene_order_Z0 <- clu_sep(CZ0_best)
cell_order_X <- clu_sep(CX_best)
cell_order_Y <- clu_sep(CY_best)
colsep <- line_sep(CZ_best);
colsep0 <- line_sep(CZ0_best);
rowsep_x <- line_sep(CX_best);
rowsep_y <- line_sep(CY_best);


## heatmap of clustering results for source data, i.e. scRNA-seq data
rowsepnum_x= clu_num(CX_best);
raw_data_x = strong_signal(X_raw[cell_order_X,gene_order_Z],CX_best,rowsepnum_x,15)
X = raw_data_x;
heatmap_fun(X, scaleyellowred, colsep, rowsep_x)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/hm_S.png "Source data")

```R
## heatmap of clustering results for the linked part of target data, i.e. scMethylation-seq data
gene_order <- clu_sep(CZ_best)
N_ave = 15;
rowsepnum_y=clu_num(CY_best);
raw_data_y = strong_signal(Y_raw_link[cell_order_Y,gene_order],CY_best,rowsepnum_y,N_ave)
heatmap_cent(raw_data_y,CZ_best)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/hm_T.png "Linked target data")

```R
## heatmap of clustering results for the unlinked part of target data
gene_order0 <- clu_sep(CZ0_best)
raw_data_y0 = strong_signal(Y_raw_unlink[cell_order_Y,gene_order0],CY_best,rowsepnum_y,N_ave)
heatmap_cent(raw_data_y0,CZ0_best)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/hm_U.png "Unlinked target data")
