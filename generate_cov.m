function covs = generate_cov(p,k_init,range)
%n number of changes, p number of variables, k number of edges changed each
%time, k_init the number at beginning.
%Kristjan Greenewald
%
%covs = cell(1,n);
%initialize
covs = .25*eye(p);
covs = add_edges(covs,k_init,range);

% for i = 2:n
%     cov_inter = remove_edges(covs{i-1},k);
%     covs{i} = add_edges(cov_inter,k,range);
% end



end

function cov = add_edges(cov,k,range)
    p = size(cov,1);
    [i,j] = rand_edges(p,k,cov == 0);
    %Weight
    a = rand*diff(range) + range(1);
    %Create edge.
    for k = 1:length(i)
        cov(i(k),i(k)) = cov(i(k),i(k)) + a;
        cov(j(k),j(k)) = cov(j(k),j(k)) + a;
        cov(i(k),j(k)) = cov(i(k),j(k)) - a;
        cov(j(k),i(k)) = cov(j(k),i(k)) - a;
    end
end

function cov = remove_edges(cov,k)
    p = size(cov,1);
    %Delete random edges.
    [i,j] = rand_edges(p,k,cov ~= 0);
    for k = 1:length(i)
        %Remove edge.
        a = -cov(i(k),j(k));
        cov(i(k),i(k)) = cov(i(k),i(k)) - a;
        cov(j(k),j(k)) = cov(j(k),j(k)) - a;
        cov(i(k),j(k)) = cov(i(k),j(k)) + a;
        cov(j(k),i(k)) = cov(j(k),i(k)) + a;
    end
end

function [i,j] = rand_edges(p,k,mask)
    mat = (ones(p) - tril(ones(p))).*mask;
    inds = randperm(nnz(mat),k);
    
    vc = mat(:);
    ix = find(vc == 1);
    vec_edges = zeros(p^2,1);
    vec_edges(ix(inds)) = 1;
    mat_edges = reshape(vec_edges,p,p);
    [i,j] = find(mat_edges);
end



