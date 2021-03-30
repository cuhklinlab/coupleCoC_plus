function [Cx, Cz, cluster_p, obj] = CoC(X,nrowcluster1,ncolcluster,iter)
%X=Y_unlink;nrowcluster1=2;ncolcluster=3;iter=15;

obj = zeros(iter,1);
p = X/sum(sum(X));%p = (p-min(p(:)))/(max(p(:))-min(p(:)));

% initialize clustering functions Cx, Cy, Cz;
rng(1);
%[Cx c] = kmeans(p, nrowcluster1);
%[Cy c] = kmeans(q, nrowcluster);
%[Cz c] = kmeans(p', ncolcluster); % p or q??

Cx = randsample(nrowcluster1, size(p,1),true); 
Cz = randsample(ncolcluster, size(p,2),true); 

% initialize tilde_p(X,Z) and tilde_q(Y,Z)
[tilde_p, cluster_p] = updateTildep_pc(p, Cx, Cz);

%tic
for t = 1:iter
  % update Cx(X) based on p, tilde_p
  Cx = updateRowClustering_cocpc(p, tilde_p, Cx);

  % update Cz(Z) based on p, q, tilde_p, tilde_q
  [Cz, dist] = updateColClustering_cocpc(p, tilde_p, Cz);

  % update tilde_p and tilde_q based on p(X,Z), Cx, Cz and q(Y,Z), Cy, Cz
  [tilde_p, cluster_p] = updateTildep_pc(p, Cx, Cz);
  
  % record values of objective function in each update
  obj(t) = dist;
end
%toc