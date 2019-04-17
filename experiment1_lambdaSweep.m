clear variables
%Script to run experiments with varying regularization parameter values. Compares L1, SCAD,
%and MCP regularization
%
%Author: K Greenewald
%Date: 11/6/16
%

K = 3; %Want to vary%K=3 FOR PLOTS
aspect_ratio = ones(1,K);
type = 'ER'; %'spiked'; %'ER' %Type of data

mag_ratio = ones(1,K);
nMC = 5; %Number of random trials to run to average plot curves.
N = 1; %*10;
%Regularization
types = {'L1','SCAD','MCP'};

a = 20;

if strcmp(type,'spiked')
    a = 2*3;
end

pscales = 32; 
lambda_rats = 1;

lambda_scale = linspace(.05,1.5,15);

tol = 1e-6;

mcc = zeros(length(pscales),length(lambda_scale),nMC);
frb = mcc;
ell2 = mcc;
%Sweep over p.
for ip = 1:length(pscales)
    
    ps = pscales(ip)*aspect_ratio;
    p = prod(ps);
    
    %Generate random ER graphs.
    if strcmp(type,'ER')
        range = [.1 .3]+.3;
        for k = 1:K
            k_init = floor(ps(k)/2*2);
            if 0 && k == 2
                Psi0{k} = mag_ratio(k)*generate_cov_grid(ps(k),k_init/2,range);
            else
                Psi0{k} = mag_ratio(k)*generate_cov(ps(k),k_init,range);
            end
        end
    elseif strcmp(type,'spiked');
        for k = 1:K
            theta = 0.5;
            Psi0{k} = mag_ratio(k)*((1-theta)*eye(ps(k)) + theta*blkdiag(ones(ps(k)/2),zeros(ps(k)/2)));
        end
    end
    
    %Get eigs
    eigz = 0;
    for k = 1:K
        [U{k},D] = eig(Psi0{k});
        if min(diag(D)) <= 0
            error('Non PSD Sigma');
        end
        eigz = eigz + kron(ones(prod(ps(1:k-1)),1),kron(diag(D),ones(prod(ps(k+1:end)),1)));
    end
    max(eigz)
    for m = 1:nMC
        S = cell(1,K);
        %Generate data
        for k = 1:K
            S{k} = 0;
        end
        %% Generate data using eigdecomp
        v = randn(prod(ps),N)./(sqrt(eigz)*ones(1,N));
        for k = 1:K
            for j = 1:prod(ps(1:k-1))
                for i = 1:prod(ps(k+1:end))
                    
                    ix = prod(ps(k:end))*(j-1) + (i + (0:ps(k)-1)*prod(ps(k+1:end)));
                   
                    v(ix,:) = U{k}*v(ix,:);
                end
            end
            
        end
        
        %% Generate Gram matrices
        for n = 1:N
            
            
            X = permute(reshape(v(:,n),ps(end:-1:1)),K:-1:1);
   
            
            %Produce sample covariances
            for k = 1:K
                XX = permute(X,[1:k-1,k+1:K, k]);
                mat = reshape(XX,prod(ps)/ps(k),ps(k));
                S{k} = S{k} + (mat'*mat)/(prod(ps)/ps(k));
            end
            
        end
        for k = 1:K
            S{k} = S{k}/N;
            
        end
        
        %Cycle over lambdas
       
        for i = 1:3 %If doing sweeps over type of regularization
        lambda_rats = lambda_rats(1)*[1,1,1];
        type = types{i};
        for j = 1:length(lambda_scale)
            lamVec = ones(1,K);
            lamVec(2) = lambda_rats(i);
            lambdas = lambda_scale(j)*sqrt(log(ps)./(p./ps)).*lamVec; 
            %Compute estimate
            tic;
            maxiter = 50;
            
            [ PsiH,~ ] = teralasso( S,ps,type,a, tol ,lambdas,maxiter);
            toc;
            %Evaluate performance
            %% MCC, not normalized
            %Compute errors.
            
            intrsect = 0;
            dsjnt = 0;
            tn = 0;
            fn = 0;
            for k = 1:K
                %errs(K,pp) = errs(K,pp) + nnz((Psi0{k} ~= 0) ~= (Psi{k} ~=0));
                
                edges_est = abs(PsiH{k} - triu(PsiH{k})) >= .1/100;% 0;
                edges = (Psi0{k} - triu(Psi0{k})) ~= 0;
                
                intrsect = intrsect + nnz(edges_est & edges);%tp
                dsjnt = dsjnt + nnz(edges_est & ~edges);%fp
                tn = tn + nnz(~edges_est & ~edges);
                fn = fn + nnz(edges & ~edges_est);
                
                
            end
             
            mcc(i,j,m) = (intrsect*tn-fn*dsjnt)/sqrt((intrsect+dsjnt)*(intrsect+fn)*(tn+dsjnt)*(tn+fn));
            %%
            if isnan(mcc(i,j,m))
                hi = 0;
            end
            
            %FrobInv, normalized
            ergs = zeros(1,p);
            ergs0 = zeros(1,p);
            for k = 1:K
                ergs = ergs + kron(ones(1,prod(ps(1:k-1))),kron(eig(PsiH{k}-Psi0{k})',ones(1,prod(ps(k+1:end)))));
                ergs0 = ergs0 + kron(ones(1,prod(ps(1:k-1))),kron(eig(Psi0{k})',ones(1,prod(ps(k+1:end)))));
            end
            frb(i,j,m) = sqrt(sum(ergs.^2))/sqrt(sum(ergs0.^2));%sqrt(p);
            %L2, not normalized
            ell = 0;
            ell0 = 0;
            for k = 1:K
                ell = ell + max(eig(PsiH{k} - Psi0{k}));
                ell0 = ell0 + max(eig(Psi0{k}));
            end
            ell2(i,j,m) = ell/ell0;
        end
        end
    end
end

%Compute means
MCC = mean(mcc,3);
FRB = mean(frb,3);
ELL2 = mean(ell2,3);

%Save
%save('K2RhoCurves.mat','pscales','aspect_ratio','MCC','FRB','ELL2','lambda_scale','nMC');

%Plot everything
%Types of regularizer
    leg = {'L1', 'SCAD', 'MCP'};

figure;

subplot(1,3,1);plot(lambda_scale,MCC,'o-');xlabel('$\rho$','interpreter','latex');ylabel('MCC');legend(leg,'interpreter','latex');
subplot(1,3,2);plot(lambda_scale,FRB,'o-');xlabel('$\rho$','interpreter','latex');ylabel('Relative $\|\hat{\Omega} - \Omega_0\|_F$','interpreter','latex');
subplot(1,3,3);plot(lambda_scale,ELL2,'o-');xlabel('$\rho$','interpreter','latex');ylabel('Relative $\|\hat{\Omega} - \Omega_0\|_2$','interpreter','latex');





