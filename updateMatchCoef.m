function Coef = updateMatchCoef(cluster_p, cluster_q, ind_X, ind_Y, gamma)
%normalize the cluster distribution
cluster_p0 = (cluster_p./sum(cluster_p,2))/sum(cluster_p,'all');
cluster_q0 = (cluster_q./sum(cluster_q,2))/sum(cluster_q,'all');
Coef = gamma * JSDiv_end(cluster_p0(ind_X,:),cluster_q0(ind_Y,:));
clearvars -except cluster_p0 cluster_q0 Coef