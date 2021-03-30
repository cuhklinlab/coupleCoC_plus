
function [ind_X, ind_Y, dis]=swap_label(cluster_p,cluster_q,k)
%normalize the cluster distribution
cluster_p0 = (cluster_p./sum(cluster_p,2))/sum(cluster_p,'all');
cluster_q0 = (cluster_q./sum(cluster_q,2))/sum(cluster_q,'all');

N_p=size(cluster_p0,1);
V_p = nchoosek(1:N_p,k);

N_q = size(cluster_q0,1);
V_q0 = nchoosek(1:N_q,k);

V_q = perms(V_q0(1,:));

for s=1:(size(V_q0,1)-1)
V_q = vertcat(V_q,perms(V_q0(s+1,:)));
end

Dis = zeros(size(V_p,1),size(V_q,1));
for i = 1:size(V_p,1)
    for j = 1:size(V_q,1)
        Dis(i,j) = JSDiv_end(cluster_p0(V_p(i,:),:),cluster_q0(V_q(j,:),:));
    end
end

[ind_x,ind_y]=find(Dis==min(min(Dis)),1);

ind_X = V_p(ind_x,:);
ind_Y = V_q(ind_y,:);
dis = min(min(Dis));
clearvars -except ind_X ind_Y dis