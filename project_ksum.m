function [ Uinvs,SigInvs,Xinvs ] = project_ksum( Xs,ps )
%PROJECT_KSUM Projects the inverse of a Kronecker sum matrix (stored as factor matrices Xs, where Xs{k} is of size ps(k) x ps(k)) onto a
%Kronecker sum of dimensions also given by the vector ps.
%   Detailed explanation goes here

K = length(Xs);
if K > 1 % If K = 1 no projection required
diagmat = zeros(ps(end:-1:1));
for i = 1:K
    [U,s] = eig(Xs{i});
    Uinvs{i} = U;
    Sigs{i} = diag(s);
    less = prod(ps(1:i-1));
    more = prod(ps(i+1:end));
    
    diagmat = diagmat + reshape(kron(ones(less,1),kron(Sigs{i},ones(more,1))),ps(end:-1:1));
    
end
diagmat = permute(diagmat,K:-1:1);

for i = 1:K
    pK = ps;
    pK(i) = 1;
    SigInvs{i} = zeros(ps(i),1);
    mtt = reshape(permute(1./diagmat,[(1:i-1) (i+1:K) i]),prod(pK),ps(i));
    %for j = 1:prod(pK)
        
    %    SigInvs{i} = SigInvs{i} + 1./(Sigs{i}*0 + mtt(:,j));
    %end
    SigInvs{i} = mean(mtt,1)';
    %SigInvs{i} = SigInvs{i}/prod(pK);
    if i == 1
        trA = mean(SigInvs{i});
    end
    %if i ~= 1
        SigInvs{i} = SigInvs{i} - trA*(K-1)/K;
    %end
end

for i = 1:K
    Xinvs{i} = Uinvs{i}*diag(SigInvs{i})*Uinvs{i}';


end
else
    %[U,s,V] = svd(Xs{1});
    [U,s] = eig(Xs{1});
    Uinvs{1} = U;
    SigInvs{1} = 1./diag(s);
    Xinvs{1} = inv(Xs{1});
end
end

