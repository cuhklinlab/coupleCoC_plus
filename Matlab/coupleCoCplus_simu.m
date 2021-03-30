clear
clc
close all

%input data
load('D:/HK_google/coupleCoC+_2021/data/simulated data/rna_simu_sy.mat');
load('D:/HK_google/coupleCoC+_2021/data/simulated data/atac_simu_sy.mat');
load('D:/HK_google/coupleCoC+_2021/data/simulated data/rna_clu_simu_sy.mat');
load('D:/HK_google/coupleCoC+_2021/data/simulated data/atac_clu_simu_sy.mat');
load('D:/HK_google/coupleCoC+_2021/data/simulated data/rna_simu_auxi_syv12.mat');
%tunining parameters
p=200;maxiter=30;
Eval_x = zeros(16,maxiter);Eval_y = zeros(16,maxiter);
%iter = 15; nrowcluster1 = 2;nrowcluster2 = 2;ncolcluster1 = 3; ncolcluster2=3;


for i = 13:30
tic
%i = 1;
ind = ((i-1)*p+1):(i*p); 
X = atac_simu(:,ind); Cx_truth = atac_clu_simu(:,i); % auxiliary data
Y_link = rna_simu(:,ind); Cy_truth = rna_clu_simu(:,i);  % target data
Y_unlink = log(rna_simu_auxi(:,ind)+1);

[X, Y_link, Cx_truth, Cy_truth] = removal_rowcol(X, Y_link, Cx_truth, Cy_truth);
yczero = sum(Y_unlink,1)==0;
Y_unlink(:,yczero)=[];

Y_full = [Y_link,Y_unlink];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%elasticC3 method
%lambda=2;beta=1;
%tic
[Cx, Cy, ~, ~, ~, ~, ~, ~,~] = coupleCoC_plus(X,Y_link,Y_unlink,2,2,3,5,15,2,1,1,2);
[~, ~, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy);
Eval_x(1:4,i) = Eval_tab{:,:}(:,1);Eval_y(1:4,i) = Eval_tab{:,:}(:,2);
%toc
%Eval_y(3,:)

%%STC method
%STC
[Cx1, Cy1, Cz1, cluster_p1, cluster_q1, obj1] = coupleCoC(X,Y_link,2,2,3,15,2);
[TAB_X1, TAB_Y1, Eval_tab1] = clu_eval(Cx_truth, Cy_truth, Cx1, Cy1);
Eval_x(5:8,i) = Eval_tab1{:,:}(:,1);Eval_y(5:8,i) = Eval_tab1{:,:}(:,2);

%%CoC_full for Target data
[Cy2, Cz2, cluster_p2, obj2] = CoC(Y_full,2,3,15);
[~, TAB_Y2, Eval_tab2] = clu_eval(Cy_truth, Cy_truth, Cy2, Cy2);
Eval_x(9:12,i) = Eval_tab2{:,:}(:,1);Eval_y(9:12,i) = Eval_tab2{:,:}(:,2);

%%k-means
[idx,~] = kmeans(X,2,'MaxIter',10000,'Replicates',15);
[idy,~] = kmeans(Y_full,2,'MaxIter',10000,'Replicates',15);
[TAB_X5, TAB_Y5, Eval_tab5] = clu_eval(Cx_truth, Cy_truth, idx, idy);
Eval_x(13:16,i) = Eval_tab5{:,:}(:,1);Eval_y(13:16,i) = Eval_tab5{:,:}(:,2);
toc
end

save('D:/HK_google/coupleCoC+_2021/data/simulated data/Eval_x_syv12','Eval_x');
save('D:/HK_google/coupleCoC+_2021/data/simulated data/Eval_y_syv12','Eval_y');

Eval_y(3:4,1:30)

load('data/simulated data/Eval_y_s1')
Az = [1:21,23:30];
mean(Eval_x(:,1:30),2);
mean(Eval_y(:,1:12),2);


[Eval_x(:,i),Eval_y(:,i)]

%%large datasets
clear
clc
close all

%input data
load('C:/Users/zengpengcheng/HK_google/elasticC3_Oct2020/data/simulated data/rna_simu_s11.mat');
load('C:/Users/zengpengcheng/HK_google/elasticC3_Oct2020/data/simulated data/atac_simu_s11.mat');
load('C:/Users/zengpengcheng/HK_google/elasticC3_Oct2020/data/simulated data/rna_clu_simu_s11.mat');
load('C:/Users/zengpengcheng/HK_google/elasticC3_Oct2020/data/simulated data/atac_clu_simu_s11.mat');
load('C:/Users/zengpengcheng/HK_google/elasticC3_Oct2020/data/simulated data/rna_simu_auxi_s11.mat');


i=1;p=1000;
ind = ((i-1)*p+1):(i*p); 
X = atac_simu(:,ind); Cx_truth = atac_clu_simu(:,i); % auxiliary data
Y_link = rna_simu(:,ind); Cy_truth = rna_clu_simu(:,i);  % target data
Y_unlink = rna_simu_auxi(:,ind);

[X, Y_link, Cx_truth, Cy_truth] = removal_rowcol(X, Y_link, Cx_truth, Cy_truth);
yczero = find(sum(Y_unlink,1)==0);
Y_unlink(:,yczero)=[];

Y_full = [Y_link,Y_unlink];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%elasticC3 method
%lambda=2;beta=1;
tic
[Cx, Cy, Cz, Cz0, cluster_p, cluster_q, cluster_q0, obj,matm] = eC3_v3(X,Y_link,Y_unlink,2,2,3,3,1,2.5,1,1,2);
[TAB_X, TAB_Y, Eval_tab] = clu_eval(Cx_truth, Cy_truth, Cx, Cy)
toc