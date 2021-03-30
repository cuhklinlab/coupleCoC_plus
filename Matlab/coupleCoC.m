function [Cx, Cy, Cz, cluster_p, cluster_q, obj] = coupleCoC_simplified(X,Y,nrowcluster1,nrowcluster2,ncolcluster,iter,lambda)
%Y = Y_link; nrowcluster1 = 2; nrowcluster2 = 3; ncolcluster = 8; iter =15; lambda = 2;

obj = zeros(iter,1);
p = X/sum(sum(X));p = (p-min(p(:)))/(max(p(:))-min(p(:)));
q = Y/sum(sum(Y));q = (q-min(q(:)))/(max(q(:))-min(q(:)));

% initialize clustering functions Cx, Cy, Cz;
rng(1);
%[Cx c] = kmeans(p, nrowcluster);
%[Cy c] = kmeans(q, nrowcluster);
%[Cz c] = kmeans(p', ncolcluster); % p or q??

Cx = randsample(nrowcluster1, size(p,1),true); 
Cy = randsample(nrowcluster2, size(q,1),true); 
Cz = randsample(ncolcluster, size(p,2),true); 

% initialize tilde_p(X,Z) and tilde_q(Y,Z)
[tilde_p, cluster_p] = updateTildep_pc(p, Cx, Cz);
[tilde_q, cluster_q] = updateTildep_pc(q, Cy, Cz);


%tic
for t = 1:iter
  % update Cx(X) based on p, tilde_p
  Cx = updateRowClustering_pc(p, tilde_p, Cx);

  % update Cy(Y) based on q, tilde_q
  Cy = updateRowClustering_pc(q, tilde_q, Cy);

  % update Cz(Z) based on p, q, tilde_p, tilde_q
  [Cz, dist] = updateColClustering_pc(p, q, tilde_p, tilde_q, Cz, lambda);

  % update tilde_p and tilde_q based on p(X,Z), Cx, Cz and q(Y,Z), Cy, Cz
  [tilde_p, cluster_p] = updateTildep_pc(p, Cx, Cz);
  [tilde_q, cluster_q] = updateTildep_pc(q, Cy, Cz);
  
  % record values of objective function in each update
  obj(t) = dist;
end
%toc