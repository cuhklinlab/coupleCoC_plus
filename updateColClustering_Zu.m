function [Cz, dist] = updateColClustering_Zu(p, tilde_p, Cz)
%p=q0;tilde_p = tilde_q0;Cz= Cx;

% compute p(X|z)
pXz = p./repmat(sum(p,1), size(p,1), 1);

% compute tilde_p(X|tilde_z)
tilde_pXz = tilde_p./repmat(sum(tilde_p,1), size(tilde_p,1), 1);

for i = 1:size(tilde_p,1)
  tilde_pXtz(i,:) = accumarray(Cz, tilde_p(i,:)')';
end

% find Cz minimizing objective function
pz = sum(p,1); 
for zc = 1:size(tilde_pXtz,2)
  for z = 1:size(pXz,2)
    temp_X(zc, z) = pz(z) * KLDiv(pXz(:,z)', tilde_pXtz(:,zc)'); 
  end
end

temp = temp_X ;

[mindist, Cz] = min(temp);
Cz = Cz';
dist = sum(mindist);

clearvars -except Cz dist
