# A quick guide to coupleCoC+ and benchmark methods

## 1. UMAP visualization for the raw data

```R
## UMAP visualization for the datasets in example 3 in the paper
# You need to install the following R packages first
library('R.matlab')
library(umap)
library(RColorBrewer)
source("R/heatmap.R")
colors = brewer.pal(10, "Set3")
# load data
X<-readMat("data/graph data/methy_X.mat");X<-as.matrix(X[[1]]);
Y<-readMat("data/graph data/methy_Y_link.mat");Y<-as.matrix(Y[[1]]);
Cx_truth<-readMat("data/graph data/methy_Cx_truth.mat");Cx_truth<-as.matrix(Cx_truth[[1]]);
Cy_truth<-readMat("data/graph data/methy_Cy_truth.mat");Cy_truth<-as.matrix(Cy_truth[[1]]);
# plot source data
Cx_truth[which(Cx_truth==1)]= rep("L4", length(which(Cx_truth==1)))
Cx_truth[which(Cx_truth==2)]= rep("L2/3 IT", length(which(Cx_truth==2)))
colors_CX = colors[3:4]
umap_x = umap(X)
plot.umap(umap_x, labels=as.factor(Cx_truth), main = "A UMAP visualization of scRNA-seq data", colors = colors_CX, cex=0.4, cex.main=1.2, cex.legend=1.2)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/ex3_S.png "Source data")

```R
# target data
Cy_truth[which(Cy_truth==1)]= rep("L4", length(which(Cy_truth==1)))
Cy_truth[which(Cy_truth==2)]= rep("L2/3 IT", length(which(Cy_truth==2)))
colors_CY = colors[3:4]
umap_y = umap(Y)
plot.umap(umap_y, labels=as.factor(Cy_truth), main = "A UMAP visualization of sc-methylation data", colors = colors_CY, cex=0.4, cex.main=1.2, cex.legend=1.2)
```
![alt text](https://github.com/cuhklinlab/coupleCoC_plus/blob/main/images/ex3_T.png "Source data")

## 2. Implementation of coupleCoC+

```MATLAB
clear
clc
close all

%%load the processed data in example 3 in the paper; note that the folder 
%"data" include all the preprocessed datasets (ex1-ex4) used in the paper.
load('data/real data/ex3_S.mat');load('data/real data/ex3_T.mat');load('data/real data/ex3_U.mat');
load('data/real data/ex3_S_cell_label.mat');load('data/real data/ex3_T_cell_label.mat'); 

%%%%%%%%%%%%%%%%%%%%% coupleCoC+ algorithm %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% input %%%%%%%%%%%%%%%
% p, q and q0 represent data S, data T and data U, respectively;
% nrowcluster1 and nrowcluster2 represent the number of cell cluster in data S, the number of cell cluster in data T, respectively;
% ncolcluster and ncolcluster0 represent the number of feature cluster in linked features and the number of feature cluster in unlinked features, respectively;
% iter: the number of iteration used in the algorithm. The default value is set as 20;
% lambda: the hyperparameter that controls the contribution of the source data S;
% beta: the hyperparameter that controls the contribution of the unlinked features in the target data T;
% gamma: the hyperparameter that controls the contribution of cell types matching across the source data S and the target data T;
% nsub: the number of matched cell clusters across source data and targert data;
%%%%% output %%%%%%%%%%%%%%%%
% Cx and Cy represent the cell label assignments for source data and target data, respectively;
% Cz and Cz0 represent the feature label assignments for the linked features and unlinked features, respectively;
% cluster_p, cluster_q and cluster_q0 represent joint probability distributions for cell clusters and features clusters in data S, data T and data U, respectively;
% obj: a vector used to record the value of objective function in each iteration;
% matm: a matrix used to record the matched order of labels in two datasets for each iteration. For example, the row [1 2 2 1] means the first and the second kind of 
% cell types in the source dataset are matched with the the second and the first kind of cell types in the target dataset, respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setting the values of hyperparameters
nrowcluster1=2;nrowcluster2=2;ncolcluster=5;ncolcluster0=8;iter=20;
lambda=0.1;beta=0.6;gamma=1;nsub=2;

% running the coupleCoC+ algorithm
[Cx, Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0, obj, matm] = coupleCoC_plus(p,q,q0,nrowcluster1,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta,gamma,nsub);

%% presenting the results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy); %% Note that this function produces the contingency table and four metrics values, two of which, i.e. ARI and NMI are utilized in this paper
disp(matm)
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
gvalue 
%% We set the number of matched N_sub as 2, because the values of the objective function g(N_sub) for 
%% choosing N_sub are smaller when N_sub = 2 (0.138 when N_sub = 1 and 0.061 when N_sub = 2, respectively)
```
```MATLAB
>> gvalue

gvalue =

    0.1382    0.0613
```

```MATLAB
%%%%%%%%%%%%%%% Tuning the hyperparameters %%%%%%%%%%%%
%%Note that we use the package in python to calculate the CH-index value for
%%each combination of hyperparameters as follows:
%%Step 1: For each combination of hyperparameters, we calculate the clustering results by coupleCoC+ and save the cell label
%%        assignments of target data; (In practice, we can tune one single hyperparameter by fixing the remaining hyperparameters);
%%Step 2: Input the data and cell labels assignments and obtain the CH-index for each combination of hyperparameters;
%%Step 3: Choose the combination that has the highest CH-index value.

% load the pacakge
from sklearn import metrics
from sklearn.metrics import pairwise_distances
import numpy as np

% calculate the CH-index
chindex = metrics.calinski_harabasz_score(data, labels)

```

## 3. Heatmaps of clustering results by coupleCoC+

```R
## import libraries
source("R/heatmap.R")
library('gplots')
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(50)
library('R.matlab')

## load the raw data and clustering results by coupleCoC+
X0<-readMat("data/graph data/methy_X.mat");X0<-as.matrix(X0[[1]]);
Y_link<-readMat("data/graph data/methy_Y_link.mat");Y_link<-as.matrix(Y_link[[1]]);
Y_unlink<-readMat("data/graph data/methy_Y_unlink.mat");Y_unlink<-as.matrix(Y_unlink[[1]]);

CX_best<-readMat("data/graph data/methy_Cx.mat");CX_best<-as.vector(CX_best[[1]])
CY_best<-readMat("data/graph data/methy_Cy.mat");CY_best<-as.vector(CY_best[[1]])
CZ_best<-readMat("data/graph data/methy_Cz.mat");CZ_best<-as.vector(CZ_best[[1]])
CZ0_best<-readMat("data/graph data/methy_Cz0.mat");CZ0_best<-as.vector(CZ0_best[[1]])

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

## 4. One example for using R-based benchmark methods in the paper
```R
## load the data
load("data/real data/dataS_ex3.Rdata")
load("data/real data/dataT_ex3.Rdata")
rna_data_spar <- rna_data_ex1 #source data
met_data_spar <- atac_data_ex1 #target data including linked part and unlinked part
colnames(rna_data_spar) <- paste("rna",1:ncol(rna_data_spar),sep='_')
colnames(met_data_spar) <- paste("met",1:ncol(met_data_spar),sep='_')
```

### LIGER
```R
# LIGER: rna+met
library(liger)

## implementation of liger
rna_met_create <- createLiger(list(rna_data=rna_data_spar,met_data=met_data_spar))
rna_met_normal <- normalize(rna_met_create)
rna_met_sele_gene <- selectGenes(rna_met_normal,datasets.use=c(1),var.thresh = 0.1,combine = "union") #only use rna to select features
rna_met_scale <- scaleNotCenter(rna_met_sele_gene)
rna_met_scale@scale.data[[2]] <- max(rna_met_scale@raw.data[[2]]) - t(as.matrix(rna_met_scale@raw.data[[2]][rna_met_scale@var.genes,])) # reverse the direction of the methylation data
rna_met_nmf = optimizeALS(rna_met_scale,k=40,remove.missing=T) #non-negative matrix factorization
rna_met_nmf_quantile2 = quantile_norm(rna_met_nmf,do.center=T)
rna_met_nmf_quantile2_louvain <- louvainCluster(rna_met_nmf_quantile2, resolution = 0.05)
result_liger <- rna_met_nmf_quantile2_louvain@clusters # this one is used in paper
```

### scACE
```R
library(mixtools)
library(label.switching)
library(scACE)
result_scACE <- getClusterGibbs(data_acc=t(rna_data_spar), data_exp=t(met_data_spar), overlap_seq_acc=1:998, overlap_seq_exp=1:1000, nCluster=3, niter=1000)
```

### Seurat
```R
library(Seurat)
##scRNA-seq data
rna_create <- CreateSeuratObject(counts = rna_data_spar)#, min.cells = 3, min.features = 200)
rna_create$tech <- "rna"
rna_norm <- NormalizeData(rna_create, normalization.method = "LogNormalize", scale.factor = 10000)
rna_feature_sele <- FindVariableFeatures(rna_norm, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(rna_feature_sele)
rna_scale <- ScaleData(rna_feature_sele, features = all.genes)
rna_pca <- RunPCA(rna_scale, features = VariableFeatures(object = rna_scale))
rna_pca_neighbor <- FindNeighbors(rna_pca, dims = 1:10)
rna_pca_neighbor_cluster <- FindClusters(rna_pca_neighbor, resolution = 0.5)
rna_seurat_label <- Idents(rna_pca_neighbor_cluster)
Cx = as.integer(rna_seurat_label);
rna_umap <- RunUMAP(rna_pca_neighbor_cluster, dims = 1:10)

##sc-Methylation data
rna_create2 <- CreateSeuratObject(counts = met_data_spar)#, min.cells = 3, min.features = 200)
rna_create2$tech <- "met"
rna_norm2 <- NormalizeData(rna_create2, normalization.method = "LogNormalize", scale.factor = 10000)
# disp / mvp / vst
rna_feature_sele2 <- FindVariableFeatures(rna_norm2, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(rna_feature_sele), 10)
all.genes2 <- rownames(rna_feature_sele2)
rna_scale2 <- ScaleData(rna_feature_sele2, features = all.genes2)
rna_pca2 <- RunPCA(rna_scale2, features = VariableFeatures(object = rna_scale2))
rna_pca_neighbor2 <- FindNeighbors(rna_pca2, dims = 1:10)
rna_pca_neighbor_cluster2 <- FindClusters(rna_pca_neighbor2, resolution = 0.5)
rna_seurat_label2 <- Idents(rna_pca_neighbor_cluster2)
Cy = as.integer(rna_seurat_label2);

## Joint analysis
### The [[]] operator can add columns to object metadata
# add true cell label to seurat object
rna_umap[["celltype"]] <- as.numeric(Cx_truth)
# identify anchors between the batch 1 scRNA dataset and the batch 2 scRNA-seq dataset and use these anchors to transfer the celltype labels.
transfer.anchors <- FindTransferAnchors(reference = rna_umap, query = rna_umap2, features = VariableFeatures(object = rna_umap), reference.assay = "RNA", query.assay = "MET", reduction = "cca")
# To transfer the cluster ids, we provide a vector of previously annotated cell type labels for the RNA to the refdata parameter. The output will contain a matrix with predictions and confidence scores for each ATAC-seq cell.
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = as.factor(rna_umap$celltype), weight.reduction = "cca")
rna_umap2 <- AddMetaData(rna_umap2, metadata = celltype.predictions)

# view the predicted cell types on a UMAP representation of the batch2 data
rna_joint_seurat_label <- celltype.predictions$predicted.id

```

### SAME-clustering
```R
library("SAMEclustering")
cluster.result_x <- individual_clustering(inputTags = rna_data_spar, datatype = "count", 
                                          percent_dropout = 10, SC3 = TRUE, CIDR = TRUE, nPC.cidr = NULL, Seurat = FALSE, nPC.seurat = NULL, 
                                          resolution = 0.9, tSNE = TRUE, dimensions = 2, perplexity = 30, SIMLR = TRUE, diverse = TRUE, SEED = 123)
cluster.ensemble_x <- SAMEclustering(Y = t(cluster.result_x), rep = 3, SEED = 123);
CX = cluster.ensemble_x$BICcluster;
```

### SIMLR
```R
library("SIMLR")
res= SIMLR_Large_Scale(X = rna_data_spar, c = 2, kk= 10) # we need to tune the parameter kk
result_simlr = res$y$cluster
```

### SC3
```R
library('rlang')
library(SingleCellExperiment)
library(SC3)

data_acc <- rna_data_spar;
sce <- SingleCellExperiment(
  assays = list(
    counts = data_acc,
    logcounts = log2(data_acc + 1)
  ), 
  colData = true_acc
)

# define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
# remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
# define spike-ins
isSpike(sce, "ERCC") <- grepl("ERCC", rowData(sce)$feature_symbol)
sce <- sc3(sce, ks = 2, gene_filter = FALSE, biology = TRUE)
result_SC3 <- as.vector(sce$sc3_2_clusters))
```

### SHARP
```R
library(cluster)
library("doParallel")
library("Matrix")
library('clusterCrit')
library('SHARP')
N = 2;
res = SHARP(rna_data_spar , logflag = FALSE, N.cluster = N, enpN.cluster=N,indN.cluster=N)
result_sharp = res$pred_clusters
```

## 5. One example for simulation study by coupleCoC+, coupleCoC and CoC
```MATLAB
%% input data for the setting 1 in simulation study
load('data/simulation data/rna_simu_s1.mat');
load('data/simulation data/atac_simu_s1.mat');
load('data/simulation data/rna_clu_simu_s1.mat');
load('data/simulation data/atac_clu_simu_s1.mat');
load('data/simulation data/rna_simu_auxi_s1v1.mat');

p=1000;maxiter=30;
Eval_x = zeros(16,maxiter);Eval_y = zeros(16,maxiter);

for i = 1:30
ind = ((i-1)*p+1):(i*p); 
X = atac_simu(:,ind); Cx_truth = atac_clu_simu(:,i); % auxiliary data
Y_link = rna_simu(:,ind); Cy_truth = rna_clu_simu(:,i);  % target data
Y_unlink = log(rna_simu_auxi(:,ind)+1);

% remove rows and columns that have zero sums
[X, Y_link, Cx_truth, Cy_truth] = removal_rowcol(X, Y_link, Cx_truth, Cy_truth);
yczero = sum(Y_unlink,1)==0;Y_unlink(:,yczero)=[];
Y_full = [Y_link,Y_unlink];

% coupleCoC+
[Cx, Cy, ~, ~, ~, ~, ~, ~,~] = coupleCoC_plus(X,Y_link,Y_unlink,2,2,3,3,15,2,1,1,2);
[~, ~, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy);
Eval_x(1:4,i) = Eval_tab{:,:}(:,1);Eval_y(1:4,i) = Eval_tab{:,:}(:,2);

% coupleCoC
[Cx1, Cy1, Cz1, cluster_p1, cluster_q1, obj1] = coupleCoC(X,Y_link,2,2,3,15,2);
[TAB_X1, TAB_Y1, Eval_tab1] = clu_eval(Cx_truth, Cy_truth, Cx1, Cy1);
Eval_x(5:8,i) = Eval_tab1{:,:}(:,1);Eval_y(5:8,i) = Eval_tab1{:,:}(:,2);

% CoC
[Cy2, Cz2, cluster_p2, obj2] = CoC(Y_full,2,3,15);
[~, TAB_Y2, Eval_tab2] = clu_eval(Cy_truth, Cy_truth, Cy2, Cy2);
Eval_x(9:12,i) = Eval_tab2{:,:}(:,1);Eval_y(9:12,i) = Eval_tab2{:,:}(:,2);
end
```
```MATLAB
>> [mean(Eval_x,2),mean(Eval_y,2)]

ans =

    1.0000    0.9827
    1.0000    0.9662
    1.0000    0.9324
    1.0000    0.9001
    0.9983    0.9733
    0.9968    0.9543
    0.9936    0.9087
    0.9918    0.8739
    0.9500    0.9500
    0.9049    0.9049
    0.8098    0.8098
    0.7659    0.7659
```


