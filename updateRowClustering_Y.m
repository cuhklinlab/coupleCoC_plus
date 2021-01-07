function [Cy] = updateRowClustering_Y(q, tilde_q, q0, tilde_q0, Cy, beta, Coef)

% compute q(Z|y) and q0(Z|y)
qZy = q./repmat(sum(q,2), 1, size(q,2));
q0Zy = q0./repmat(sum(q0,2), 1, size(q0,2));

% compute tilde_q(Z|tilde_y) and tilde_q0(Z|tilde_y)
tilde_qZy = tilde_q./repmat(sum(tilde_q,2), 1, size(tilde_q,2));
tilde_q0Zy = tilde_q0./repmat(sum(tilde_q0,2), 1, size(tilde_q0,2));

for i = 1:size(tilde_q,2)
  tilde_qZty(:,i) = accumarray(Cy, tilde_qZy(:,i));
end

for i = 1:size(tilde_q0,2)
  tilde_q0Zty(:,i) = accumarray(Cy, tilde_q0Zy(:,i));
end

% find Cx minimizing D(p(Z|x) || tilde_p(Z|tilde_x))

for rc = 1:size(tilde_qZty, 1)
  for r = 1:size(qZy, 1)
    temp(rc, r) = KLDiv(qZy(r,:), tilde_qZty(rc,:)) + beta * KLDiv(q0Zy(r,:), tilde_q0Zty(rc,:)) + Coef / (sum(q(r,:),2) * size(q,1));
  end
end

[mindist Cy] = min(temp);
Cy = Cy';

clearvars -except Cy
