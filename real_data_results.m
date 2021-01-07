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
