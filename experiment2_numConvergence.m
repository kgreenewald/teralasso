%Script to illustrate numerical convergence of TeraLasso algorithm
clear variables
type = 'MCP';%Type of regularization
rng(33);
opt = 1;
K = 4;
big = 1;
errNrm = 'fro';

ps = 40*ones(1,K);
N = 1;
a = 60;

tol = 1e-10;
lk = 1;
lambda = lk*ones(1,K)/(sqrt(N*prod(ps)/ps(1)));
if N==1
    lambda = lambda/3;
end
range = .1+[.1 .3];
mn = 0;
for k = 1:K
    k_init = floor(ps(k));
    Psi0{k} = generate_cov(ps(k),k_init,range);
    mn = mn + mean(diag(Psi0{k}));
end

%Find eigenvalues

eigz = 0;
for k = 1:K
    [U{k},D] = eig(Psi0{k});
    if min(diag(D)) <= 0
        error('Non PSD Sigma');
    end
    eigz = eigz + kron(ones(prod(ps(1:k-1)),1),kron(diag(D),ones(prod(ps(k+1:end)),1)));
end
for k = 1:K
    S{k} = 0;
end
%Generate data using eigdecomp
v = randn(prod(ps),N)./(sqrt(eigz)*ones(1,N));
for k = 1:K
    for j = 1:prod(ps(1:k-1))
        for i = 1:prod(ps(k+1:end))
            
            ix = prod(ps(k:end))*(j-1) + (i + (0:ps(k)-1)*prod(ps(k+1:end)));
            
            v(ix,:) = U{k}*v(ix,:);
        end
    end
    
end

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


maxiter = 200;
tic

[ Psi,Psis ] = teralasso( S,ps,type,a,tol ,lambda,maxiter);
time4 = toc;

stat_err4 = 0;
for k = 1:K
    diff = (Psi0{k} - Psi{k}) - eye(ps(k))*mean(diag((Psi0{k} - Psi{k})));
    if errNrm ~= inf
        stat_err4 = stat_err4 + prod(ps)/ps(k)*norm(diff,errNrm)^2;
    else %infinity norm
        stat_err4 = max(stat_err4, max(max(abs(diff))));
    end
end
for k = 1:K
    Psis{1}{k} = eye(ps(k));
end
for count = 1:length(Psis)
    nrm(count) = 0;
    for k = 1:K
        diff = (Psis{count}{k} - Psi{k}) - eye(ps(k))*mean(diag((Psis{count}{k} - Psi{k})));
        if errNrm ~= inf
            nrm(count) = nrm(count) + prod(ps)/ps(k)*norm(diff,errNrm)^2;
        else %infinity norm
            nrm(count) = max(nrm(count), max(max(abs(diff))));
        end
    end
end
%% PLOT
ix = find(nrm./nrm(1)< 1e-6,1,'first');
if N == 1
    ix = find(nrm./nrm(1)< 1e-4,1,'first');
end
ix = find(nrm < stat_err4,1,'first')+2;%FOR NONCONVEX
figure;hold off; semilogy(sqrt(nrm(1:ix-1)./nrm(1)),'o-b');
hold on;semilogy(1:length(nrm(1:ix-1)),sqrt(stat_err4/nrm(1))*ones(1,length(nrm(1:ix-1))),'k');
ixx = find(nrm < stat_err4,1,'first');
text(1, 1.2*sqrt(stat_err4/nrm(1)),[num2str(ixx/length(nrm)*time4,3) ' sec']);


K = 2;
ps = ps(1)^opt*ones(1,K);

lambda = lk*ones(1,K)/(sqrt(N*prod(ps)/ps(1)));


if opt == 1
    v = v(1:prod(ps),:);
end

%% Get S
for k = 1:K
    S{k} = 0;
end
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

%% Estimate
tic

[ Psi,Psis ] = teralasso( S,ps,type,a,tol ,lambda,maxiter);

time2 = toc;
stat_err2 = 0;

for k = 1:K
    if opt == 2
        mm = kron(Psi0{(k-1)*2 + 1},eye(size(Psi0{(k-1)*2 + 2})))+ kron(eye(size(Psi0{(k-1)*2 + 1})),Psi0{(k-1)*2 + 2});
    else
        mm = Psi0{k};
    end
    diff = (mm - Psi{k}) - eye(ps(k))*mean(diag((mm - Psi{k})));
    
    if errNrm ~= inf
        stat_err2 = stat_err2 + prod(ps)/ps(k)*norm(diff,errNrm)^2;
    else %infinity norm
        stat_err2 = max(stat_err2, max(max(abs(diff))));
    end
end
for k = 1:K
    Psis{1}{k} = eye(ps(k));
end
for count = 1:length(Psis)
    nrm(count) = 0;
    for k = 1:K
        diff = (Psis{count}{k} - Psi{k}) - eye(ps(k))*mean(diag((Psis{count}{k} - Psi{k})));
        
        if errNrm ~= inf
            nrm(count) = nrm(count) + prod(ps)/ps(k)*norm(diff,errNrm)^2;
        else %infinity norm
            nrm(count) = max(nrm(count), max(max(abs(diff))));
        end
    end
end
if opt == 1
    nrm = nrm*100;
end

%% Plot
ix = find(nrm./nrm(1)< 1e-4,1,'first');
if N == 1
    ix = find(nrm./nrm(1)< .2e-2,1,'first');
end
ix = find(nrm < stat_err2,1,'first')+2;
hold on;semilogy(sqrt(nrm(1:ix-1)./nrm(1)),'rx-');
hold on;semilogy(1:length(nrm(1:ix-1)),sqrt(stat_err2/nrm(1))*ones(1,length(nrm(1:ix-1))),'g');
ixx = find(nrm < stat_err2,1,'first');
text(1, 1.2*sqrt(stat_err2/nrm(1)),[num2str(ixx/length(nrm)*time2,3) ' sec']);
xlabel('Iterations');
ylabel('Norm. $\|\Omega_t - \Omega^*\|_F$','interpreter','latex')
if opt == 2
    if big == 1
        
        legend({'K = 4, p = 2.56x10^6','K = 4 stat. error',  'K = 2, p = 2.56x10^6','K = 2 stat. error'});
    else
        legend({'K = 4, p = 10^4','K = 4 stat. error', 'K = 2, p = 10^4','K = 2 stat. error'});
    end
else
    if big == 0
        legend({'K = 4, d_k = 10','K = 4 stat. error', 'K = 2, d_k = 10','K = 2 stat. error'});
    else
        
        legend({'K = 4, d_k = 40','K = 4 stat. error', 'K = 2, d_k = 40','K = 2 stat. error'});
    end
end
title(['Convergence of TeraLasso with ' type ' regularization'])