function [X, Y, Cx_truth, Cy_truth] = removal_rowcol(X, Y, Cx_truth, Cy_truth)
% find the colums and rows with zero sum in X and Y;
czero = [find(sum(X,1)==0),find(sum(Y,1)==0)];
yrzero = find(sum(Y,2)==0); xrzero = find(sum(X,2)==0);
Y(:,czero)=[];Y(yrzero,:)=[];
X(:,czero)=[];X(xrzero,:)=[];
Cx_truth(xrzero)=[];Cy_truth(yrzero)=[];