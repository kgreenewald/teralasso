function [ auc ] = auc_comp( likeW,likeS )
%AUC Summary of this function goes here
%   Detailed explanation goes here

lWSort = sort(likeW,'ascend');
for i = 1:length(likeW);
    pfa(i) = (i-1)/length(likeW);
    pd(i) = nnz(likeS <= lWSort(i))/length(likeS);
end

auc = trapz(pfa,pd);

auc = max(auc,1-auc);

end

