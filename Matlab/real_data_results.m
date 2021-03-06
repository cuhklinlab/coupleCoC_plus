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

%%results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy);
disp(matm)

%%determination of N_sub
gvalue = [];
for N_sub=1:2
[~, ~, dis] = swap_label_plus(cluster_p,cluster_q,N_sub);
gvalue(N_sub)= dis/(N_sub*log(N_sub+1));
end
gvalue










clear
clc
close all

% batch effect analysis for example 4
load('data/celltype_b1.mat');load('data/celltype_b2.mat');
load('data/data_batch1.mat');load('data/data_batch2.mat'); 

data_b1_sub = data_batch1(1:1000,98:385);%find(celltype_b1==1)
Cx_truth = celltype_b1(98:385);
data_b2_sub = data_batch2(1:1000,1:288);%find(celltype_b2==4)
Cy_truth = celltype_b2(1:288);

p = log(data_b1_sub+1).';
q = log(data_b2_sub+1).';

%%coupleCoC+
%setting hyperparameters
nrowcluster1=3;nrowcluster2=3;ncolcluster=5;ncolcluster0=5;iter=20;
lambda=2;beta=0;gamma=1;nsub=1;
[Cx, Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0, obj, matm] = coupleCoC_plus(p,q,q,nrowcluster1,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta,gamma,nsub);

%%results
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy);
disp(matm)

%
save('batch1_cluster_p','cluster_p');save('batch2_cluster_q','cluster_q');save('batch_matm','matm');
load('D:/HK_google/coupleCoC+_2021/Batcheffect/Cx.mat');
load('D:/HK_google/coupleCoC+_2021/Batcheffect/Cy.mat');

N_sub = 3;
[ind_X, ind_Y, dis] = swap_label_plus(cluster_p,cluster_q,N_sub);

%%batch mixing
Cbm_truth = [celltype_b1(98:289).',celltype_b2(97:288).'];
%Cbm_truth = [repelem(1,192),repelem(2,192)];
Cy0 = repelem(0,length(Cy));
Cy0(Cy==3) = repelem(1,length(find(Cy==3)));
Cy0(Cy==1) = repelem(3,length(find(Cy==1)));
Cy0(Cy==2) = repelem(2,length(find(Cy==2)));
Cbm = [Cx(1:192).',Cy0(97:288)];
[TAB_bm, TAB_bm, Eval_bm] = clu_eval(Cbm_truth, Cbm_truth, Cbm, Cbm);

save('D:/HK_google/coupleCoC+_2021/Batcheffect/Cx.mat','Cx');
save('D:/HK_google/coupleCoC+_2021/Batcheffect/Cy.mat','Cy');
save('D:/HK_google/coupleCoC+_2021/Batcheffect/Cz.mat','Cz');

save('D:/HK_google/coupleCoC+_2021/Batcheffect/X.mat','p');
save('D:/HK_google/coupleCoC+_2021/Batcheffect/Y.mat','q');

save('D:/HK_google/coupleCoC+_2021/Batcheffect/Cx_truth.mat','Cx_truth');
save('D:/HK_google/coupleCoC+_2021/Batcheffect/Cy_truth.mat','Cy_truth');