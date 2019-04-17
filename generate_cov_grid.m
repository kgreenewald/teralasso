function covs = generate_cov_grid(p,k_init,range)
%n number of changes, p number of variables, k number of edges changed each
%time, k_init the number at beginning.
%Makes a random grid graph. p must be a square number.
%Range should be [-.5, .5]. k_init should be small.
%Kristjan Greenewald
%

%initialize
covs = .25*eye(p);
covs = add_edges(covs,k_init,range);




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
    %Must fit a grid.
    
    %Grid graph - each node is connected to above and below, side and back.
    grid_edge = zeros(p);
    sz = [sqrt(p) sqrt(p)];
    for ii = 1:sqrt(p)
        for jj = 1:sqrt(p)
            ix_now = sub2ind(sz,ii,jj);
            if ii > 1
                ix = sub2ind(sz,ii-1,jj);
                grid_edge(ix_now,ix) = 1;
                grid_edge(ix,ix_now) = 1;
            end
            if jj > 1
                ix = sub2ind(sz,ii,jj-1);
                grid_edge(ix_now,ix) = 1;
                grid_edge(ix,ix_now) = 1;
            end
            if ii < sqrt(p)
                ix = sub2ind(sz,ii+1,jj);
                grid_edge(ix_now,ix) = 1;
                grid_edge(ix,ix_now) = 1;
            end
            if jj < sqrt(p)
                ix = sub2ind(sz,ii,jj+1);
                grid_edge(ix_now,ix) = 1;
                grid_edge(ix,ix_now) = 1;
            end
            
        end
    end
    
    
    %mat = (ones(p) - tril(ones(p))).*mask;
    mat = grid_edge.*mask;
    
    %Shuffle available edges
    inds = randperm(nnz(mat),min(k,nnz(mat)));
    
    %Indices of available edges.
    vc = mat(:);
    ix = find(vc == 1);
    %Form adjacency matrix
    vec_edges = zeros(p^2,1);
    %Apply selected edges.
    vec_edges(ix(inds)) = 1;
    mat_edges = reshape(vec_edges,p,p);
    %Pull out indices.
    [i,j] = find(mat_edges);
end



