function [ Psi,Psis ] = teralasso( S,ps,type,a,tol ,lambda,maxiter)
%TERALASSO Implementation of Tensor Graphical Lasso (TeraLasso) precision
%matrix estimator. Learns a K-order Kronecker sum model
%
%Inputs:
%   S: Length K cell array of factorwise Gram matrices (S_k in paper)
%   ps: Length K vector, ps(k) = size(S{k},1). (= d_k in paper)
%   type: Type of regularization: L1, SCAD, MCP.
%   a: Parameter for SCAD/MCP.
%   tol: Tolerance for convergence criterion.
%   lambda: Length K vector, consisting of L1 penalties for each factor.
%   maxiter: Maximum number of iterations allowed
%
%Author:Kristjan Greenewald

%Nonconvex contraints.
switch type
    case 'L1'
        mu = 1e-6;
    case 'SCAD'
        mu = 1/(a-1);
    case 'MCP'
        mu = 1/a;
end
kappa = sqrt(2/mu)*2;


K = length(ps);
p = prod(ps);
zeta0 = .1;
c = .5;%backtracking constant
for k = 1:K
    Psi{k} = 1/K*eye(ps(k));
    M{k} = 1/K*eye(ps(k));
end
count = 1;
nrm = inf;
lambda = lambda + 1e-9;%Not strictly required, for large, small sample size covariances when lambda = 0, the inverted matrix may be ill conditioned.



logdet = 0;
while max(nrm(max(1,count-3):count)) > tol && count < maxiter%? > tol
    count = count + 1;
    %count
    %Line search
    zeta = zeta0;
    for k = 1:K
        Sp = S{k} - (K-1)/K*mean(diag(S{k}))*eye(size(S{k}));
        %PsiN{k} = shrinkL1(Psi{k} - zeta*(Sp - M{k}),lambda(k)*zeta);
        PsiN{k} = shrink_regularizer(Psi{k} - zeta*(Sp - M{k}),lambda(k),zeta,type,a);
    end
    if zeta~=0
        cnt = 0;
        [result,result1,logdetN] = eval_cond(PsiN,Psi,M,S,ps,zeta,logdet,kappa);
        while ~result && cnt <= 10 %eval_cond(PsiN,Psi,S,ps) && cnt <= 10
            cnt = cnt + 1;
            zeta = zeta*c;
            for k = 1:K
                Sp = S{k} - (K-1)/K*mean(diag(S{k}))*eye(size(S{k}));
                %PsiN{k} = shrinkL1(Psi{k} - zeta*(Sp - M{k}),lambda(k)*zeta);
                PsiN{k} = shrink_regularizer(Psi{k} - zeta*(Sp - M{k}),lambda(k),zeta,type,a);
            end
            [result,result1,logdetN] = eval_cond(PsiN,Psi,M,S,ps,zeta,logdet,kappa);
        end
        
        if cnt == 11
            zeta = 0;
            for k = 1:K
                zeta = zeta + min(eig(Psi{k}));
            end
            zeta = zeta^2/2;
            for k = 1:K
                Sp = S{k} - (K-1)/K*mean(diag(S{k}))*eye(size(S{k}));
                PsiN1{k} = shrinkL1(Psi{k} - zeta*(Sp - M{k}),lambda(k)*zeta);
            end
            [result,result2,logdetN] = eval_cond(PsiN,Psi,M,S,ps,zeta,logdet,kappa);
            if result2
                PsiN = PsiN1;
            end
        end
        if ~isreal(logdetN)
            hi =0;
        end
        
        logdet = logdetN;
    end
    %Set next initial step
    [~,~,MsN] = project_ksum(PsiN,ps);%Project onto Kronecker sum subspace
    Num = 0;
    Den = 0;
    for k = 1:K
        dgI(k) = mean(diag((PsiN{k}-Psi{k})));
        Num = Num + prod(ps([1:k-1 k+1:K]))*sum(sum((PsiN{k}-Psi{k} - eye(ps(k))*dgI(k)).^2));
        dK = PsiN{k} - Psi{k};
        for kk = 1:K
            dKK = M{kk} - MsN{kk};
            if k == kk
                Den = Den + prod(ps([1:k-1 k+1:K]))*sum(sum(dK.*dKK));
            else
                Den = Den + trace(dK)*trace(dKK)*p/(ps(k)*ps(kk));
            end
        end
    end
    Num = Num + mean(dgI)^2*prod(ps);
    zeta0 = Num/Den;
    zet(count) = zeta0;
    if Den == 0
        zeta0 = 0;
    elseif abs(zeta0) > 1e2
        zeta0 = 1e2*sign(zeta0);
        
    elseif zeta0 < 0
        zeta0 = 1;
        
    end
    %Update iterate
    PsiOld = Psi;
    Psi = PsiN;
    M = MsN;
    
    %Compute convergence
    nrm(count) = 0;
    for k = 1:K
        diff = (Psi{k} - PsiOld{k}) - eye(ps(k))*mean(diag((Psi{k} - PsiOld{k})));
        nrm(count) = nrm(count) + norm(diff,'fro')^2;
    end
    
    for k = 1:K
        Psis{count}{k} = Psi{k};
    end
    
    
   
end

end

%%%%%%%%%%%

function [result,result1,logdetN] = eval_cond(PsiN,Psi,M,Ss,ps,zeta,logdet,kappa)
K = length(ps);
p = prod(ps);
mineig = 0;
ergNvec = 0;
%ergvec = 0;
for k = 1:K
    ergsN{k} = eig(PsiN{k});
    %ergNvec = ergNvec + kron(ones(prod(ps(1:k-1)),1),kron(ergsN{k},ones(prod(ps(k+1:end)),1)));
    ergNvec = ergNvec + repmat(kron(ergsN{k},ones(prod(ps(k+1:end)),1)),prod(ps(1:k-1)),1);
    %ergvec = kron(ones(prod(ps(1:k-1)),1),kron(ergs{k},ones(prod(ps(k+1:end)),1)));
    %mineig = mineig + min(ergsN);
end
mineig = min(ergNvec);
logdetN = -sum(log(ergNvec));
if mineig > 0 && max(ergNvec) <=kappa
    result1 = 1;
else
    result1 = 0;
end
if mineig < 0
    result = 0;
    return;
else
    %Norm of change
    Num = 0;
    
    for k = 1:K
        dgI(k) = mean(diag((PsiN{k}-Psi{k})));
        Num = Num + prod(ps([1:k-1 k+1:K]))*sum(sum((PsiN{k}-Psi{k} - eye(ps(k))*dgI(k)).^2));
        
    end
    Num = Num + mean(dgI)^2*prod(ps);
    
    
    %Rest
    %result = 1;
    
    fNew = logdetN;
    for k = 1:K
        fNew = fNew + sum(sum(Ss{k}.*PsiN{k}))*prod(ps)/ps(k);
    end
    %logdetN = -sum(log(ergvec));
    Q = logdet;
    
    for k = 1:K %
        Q = Q + sum(sum(Ss{k}.*Psi{k}))*prod(ps)/ps(k)  +0 + 1/(2*zeta)*Num;
        Sp = Ss{k} - (K-1)/K*mean(diag(Ss{k}))*eye(size(Ss{k}));
        for kk = 1:K
            if kk == k
                
                Q = Q + sum(sum((Sp-M{k}).*(PsiN{k} - Psi{k})))*prod(ps)/ps(k);
            else
                Q = Q - trace(M{k}-Sp)*trace((PsiN{kk} - Psi{kk}))*prod(ps)/(ps(k)*ps(kk));
            end
        end
    end
    
    if fNew <= Q
        result = 1;
    else
        result = 0;
    end
end

end