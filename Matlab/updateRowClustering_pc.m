function [Cx] = updateRowClustering_pc(p, tilde_p, Cx)

% compute p(Z|x)
pZx = p./repmat(sum(p,2), 1, size(p,2));

% compute tilde_p(Z|tilde_x)
tilde_pZx = tilde_p./repmat(sum(tilde_p,2), 1, size(tilde_p,2));

for i = 1:size(tilde_p,2)
  tilde_pZtx(:,i) = accumarray(Cx, tilde_pZx(:,i));
end

% find Cx minimizing D(p(Z|x) || tilde_p(Z|tilde_x))

for rc = 1:size(tilde_pZtx, 1)
  for r = 1:size(pZx, 1)
    temp(rc, r) = KLDiv(pZx(r,:), tilde_pZtx(rc,:));% + beta * JSDiv_pc(cluster_p, cluster_q, rc,indicator) / (sum(p(r,:),2) * size(p,1) * size(tilde_pZtx, 1));
  end
end

[mindist Cx] = min(temp);
Cx = Cx';

clearvars -except Cx
