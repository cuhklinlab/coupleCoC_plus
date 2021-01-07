function [TAB_X, TAB_Y, Eval_tab] = clu_eval(clu_X_truth, clu_Y_truth, clu_X_bes, clu_Y_bes)
sample = zeros(4,2);
sample(1,1)=Cal_Purity(clu_X_truth, clu_X_bes);
sample(1,2)=Cal_Purity(clu_Y_truth, clu_Y_bes);
[sample(3,1),sample(2,1),~,~]=RandIndex(clu_X_truth,clu_X_bes);
[sample(3,2),sample(2,2),~,~]=RandIndex(clu_Y_truth,clu_Y_bes);
sample(4,1) = Cal_NMI(clu_X_truth, clu_X_bes);
sample(4,2) = Cal_NMI(clu_Y_truth, clu_Y_bes);
rowNames = {'Purity','RI','ARI','NMI'};
colNames = {'X','Y'};
Eval_tab = array2table(sample,'RowNames',rowNames,'VariableNames',colNames);

TAB_X = crosstab(clu_X_truth,clu_X_bes);
TAB_Y = crosstab(clu_Y_truth,clu_Y_bes);