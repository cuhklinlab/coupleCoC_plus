function [Cx, Cy, Cz, Cu, cluster_p, cluster_q, cluster_q0, obj, matm] = coupleCoC_plus(p,q,q0,nrowcluster1,nrowcluster2,ncolcluster,ncolcluster0,iter,lambda,beta,gamma,nsub)
%p: preprocessed S;
%q: preprocessed T;
%q0: preprocessed U;
%matm: matching results

obj = zeros(iter,1);
matm = zeros(iter,2*nsub);

% initialize clustering functions Cx, Cy, Cz;
rng(1);
Cx = randsample(nrowcluster1, size(p,1),true); 
Cy = randsample(nrowcluster2, size(q,1),true); 
Cz = randsample(ncolcluster, size(p,2),true); 
Cu = randsample(ncolcluster0, size(q0,2),true); 

% initialize tilde_p(X,Z), tilde_q(Y,Z) and tilde_q0(Y,Z0)
[tilde_p, cluster_p] = updateTildep_plus(p, Cx, Cz);
[tilde_q, cluster_q] = updateTildep_plus(q, Cy, Cz);
[tilde_q0, cluster_q0] = updateTildep_plus(q0, Cy, Cu);
[ind_X, ind_Y]=swap_label_plus(cluster_p,cluster_q,nsub);


for t = 1:iter
  % upate the coefficient of the matching term D_{KL}(...||...)
  Coef = updateMatchCoef(cluster_p, cluster_q, ind_X, ind_Y, gamma);
  
  % update Cx(X) based on p, tilde_p, Coef;
  Cx = updateRowClustering_X(p, tilde_p, Cx, Coef);

  % update Cy(Y) based on q, tilde_q, q0, tilde_q0, Coef;
  Cy = updateRowClustering_Y(q, tilde_q, q0, tilde_q0, Cy, beta, Coef);
  
  % updata Cz0(Z0) based on q0, tilde_q0
  [Cu, dist0] = updateColClustering_Zu(q0, tilde_q0, Cu);

  %[Cz, dist] = updateColClustering_cocpc(p, tilde_p, Cz);
  % update Cz(Z) based on p, q, tilde_p, tilde_q
  [Cz, dist] = updateColClustering_Z(p, q, tilde_p, tilde_q, Cz, lambda, Coef);

  % update tilde_p and tilde_q based on p(X,Z), Cx, Cz and q(Y,Z), Cy, Cz
  [tilde_p, cluster_p] = updateTildep_plus(p, Cx, Cz);
  [tilde_q, cluster_q] = updateTildep_plus(q, Cy, Cz);
  [tilde_q0, cluster_q0] = updateTildep_plus(q0, Cy, Cu);
  
  % match the cell types across the target data and the source data;
  [ind_X, ind_Y]=swap_label_plus(cluster_p,cluster_q,nsub);
  matm(t,:) = [ind_X,ind_Y]; %record the matching history;
  
  % record values of objective function in each update
  obj(t) = dist + dist0 + Coef;
end
