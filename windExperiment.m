%Script to recreate sparse Figure 12 (Wind speed experiment). ML Kron Product results
%omitted to control runtime. For full figure 12 comment line 7.
clear variables
close all

FF_iter = 100;%Max number of iterations of the flip/flop BiGLasso algorithm.
FF_iter = 1; %Do not compute BiGLasso because of run length.


load('WindData.mat','spd','meanS');%Relevant portion of full dataset

%Do some covariances.

if 1
    xv = -89 + (104:1:124);%eastern grid
else
    xv = -89 + (90:109);%western grid
end
yv = -22 + (23:1:32);

f = length(xv)*length(yv);

%Make the matrix.

samps = [1 5:5:25]; 
tt= [5 10 15]; %Temporal widths of covariance

for ii = 1:3
    for k = 1:length(samps)
        n = samps(k);%365;
        for season = 0:1
            t = tt(ii);
            mt = zeros(n,length(xv)*length(yv)*t);
            SS{1} =0;
            SS{2} = 0;
            SS{3} = 0;
            
            ps = [t length(yv) length(xv)];
            K = length(ps);
            for i = 1:n
                
                mtt = spd{1}(xv,yv,364/2*season + i+(0:t-1)) - 1*meanS(xv,yv,364/2*season + i+(0:t-1));
                mt(i,:) = mtt(:)';
                
                for kk = 1:K
                    XX = permute(mtt,[1:kk-1,kk+1:K, kk]);
                    mat = reshape(XX,prod(ps)/ps(kk),ps(kk));
                    SS{kk} = SS{kk} + (mat'*mat)/(prod(ps)/ps(kk));
                end
            end
            
            scm = mt'*mt/n;
            
            %First without regularization
            lambda = zeros(1,3);
            tol = 1e-8;
            
            SS{1} = SS{1}/n;
            SS{2} = SS{2}/n;
            SS{3} = SS{3}/n;
            scmK = kron(SS{1},kron(SS{2},SS{3}));
            SK = SS;
            maxiter = 90;
            tic;[ Psi ] = teralasso( SS,ps,'L1',0,tol ,lambda,maxiter);toc;
            
            %Second with regularization %1/sqrt(t) on third axis comes from theoretically best regularizer with d_3=t.
            lambda = 3.5e-5/sqrt(n)*[1 1 1/sqrt(t)];
            
            tic;[ PsiL ] = teralasso( SS,ps,'L1',0,tol ,lambda,maxiter);toc;
            
            
            Sig = inv(kron(Psi{1},eye(ps(2)*ps(3))) + kron(eye(ps(1)),kron(Psi{2},eye(ps(3))))+kron(eye(ps(1)*ps(2)),Psi{3}));
            Omeg = kron(Psi{1},eye(ps(2)*ps(3))) + kron(eye(ps(1)),kron(Psi{2},eye(ps(3))))+kron(eye(ps(1)*ps(2)),Psi{3});
            OmegL = kron(PsiL{1},eye(ps(2)*ps(3))) + kron(eye(ps(1)),kron(PsiL{2},eye(ps(3))))+kron(eye(ps(1)*ps(2)),PsiL{3});
            
            
            [ A_mat, B_mat, C_mat] = FFKron3(scm,ps(1),ps(2),ps(3),0,FF_iter, tol);
           
            OmegP = kron(inv(A_mat),kron(inv(B_mat),inv(C_mat)));
            yix = 1 + (0:ps(2)*ps(3):t-1);
            xix = ones(1,prod(ps));
            xix(yix) = 0;
            xix = logical(xix);
            
            pred_coef = Sig(yix,xix)/(Sig(xix,xix));
            pred_coefKron = scmK(yix,xix)*pinv(scmK(xix,xix));
            pred_coefSCM = scm(yix,xix)*pinv(scm(xix,xix));
            maxI =30;
            initer = sum(log(eig(Omeg)))*maxI;
            initerL = sum(log(eig(OmegL)))*maxI;
            initerP = sum(log(eig(OmegP)))*maxI;
            tic;
            for j = (2:51)
                
                for i = 1:maxI-t+1
                    mtt = spd{j}(xv,yv,i+(0:t-1)) - 1*meanS(xv,yv,i+(0:t-1));
                    mtt = mtt(:);
                    
                    prd(:,i + (j-1)*maxI) = pred_coef*mtt(xix);
                    prdSCM(:,i + (j-1)*maxI) = pred_coefSCM*mtt(xix);
                    prdSCMK(:,i + (j-1)*maxI) = pred_coefKron*mtt(xix);
                    act(:,i+(j-1)*maxI) = mtt(yix);
                    
                    
                    
                end
                %Likelihoods
                
                likeW(j) = initer;
                likeS(j) = initer;
                likeWL(j) = initerL;
                likeSL(j) = initerL;
                likeWP(j) = initerP;
                likeSP(j) = initerP;
                for i=1:maxI-t%*1.5 %Days
                    %Winter
                    mtt = spd{j}(xv,yv,i+(0:t-1)) - 1*meanS(xv,yv,i+364/2*season+(0:t-1));
                    mtt = mtt(:);
                
                    
                    likeW(j) = likeW(j)-mtt'*Omeg*mtt;
                    likeWP(j) = likeWP(j)-mtt'*OmegP*mtt;
                    likeWL(j) = likeWL(j) - mtt'*OmegL*mtt;
                    %Summer
                    mtt = spd{j}(xv,yv,i+364/2+(0:t-1)) - 1*meanS(xv,yv,i+364/2*season+(0:t-1));
                    mtt = mtt(:);
                   
                    likeS(j) = likeS(j) - mtt'*Omeg*mtt;
                    likeSL(j) = likeSL(j) - mtt'*OmegL*mtt;
                    likeSP(j) = likeSP(j) - mtt'*OmegP*mtt;
                end
                
            end
            toc
            likeWC{season+1} = likeW;
            likeSC{season+1} = likeS;
            likeWLC{season+1} = likeWL;
            likeSLC{season+1} = likeSL;
            likeWPC{season+1} = likeWP;
            likeSPC{season+1} = likeSP;
        end
        
        m1 = mean([likeWC{1}(2:end) likeSC{1}(2:end)]);
        m2 = mean([likeWC{2}(2:end) likeSC{2}(2:end)]);
        mL1 = mean([likeWLC{1}(2:end) likeSLC{1}(2:end)]);
        mL2 = mean([likeWLC{2}(2:end) likeSLC{2}(2:end)]);
        mP1 = mean([likeWPC{1}(2:end) likeSPC{1}(2:end)])*0;
        mP2 = mean([likeWPC{2}(2:end) likeSPC{2}(2:end)])*0;
        
        err(k) = 1-(nnz(likeWC{0+1}(2:end)-m1 > likeWC{1+1}(2:end)-m2) + nnz(likeSC{0+1}(2:end)-m1 <likeSC{1+1}(2:end)-m2))/(2*length(likeWC{0+1}(2:end)));
        errL(k) = 1-(nnz(likeWLC{0+1}(2:end)-mL1 > likeWLC{1+1}(2:end)-mL2) + nnz(likeSLC{0+1}(2:end)-mL1 < likeSLC{1+1}(2:end)-mL2))/(2*length(likeWLC{0+1}(2:end)));
        errP(k) = 1-(nnz(likeWPC{0+1}(2:end)-mP1 > likeWPC{1+1}(2:end)-mP2) + nnz(likeSPC{0+1}(2:end)-mP1 < likeSPC{1+1}(2:end)-mP2))/(2*length(likeWPC{0+1}(2:end)));
        
        auc(k) = auc_comp( likeW(2:end),likeS(2:end) );
        aucP(k) = auc_comp( likeWP(2:end),likeSP(2:end) );
        aucL(k) = auc_comp( likeWL(2:end),likeSL(2:end) );
        
    end
    if FF_iter == 1 %Don't show FF plot
        errP = NaN*errP;
    end
    figure(4);subplot(1,3,ii);hold off; plot(samps,err,'bo-');hold on;plot(samps,errL,'ko-');hold on;plot(samps,errP,'ro-');
    xlabel('n');
    ylabel('Error Rate');
    legend('Proposed ML Kron. Sum','Proposed TeraLasso', 'ML Kron. Product (FF)');
    title(['T = ' num2str(t)]);
    axis([1 25 0 .5]);
end
