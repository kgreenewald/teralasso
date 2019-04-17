function [ A_mat, B_mat, C_mat] = FFKron3(SCM,p,f,g,n,niter, tol)
%FF Function implements the flip-flop algorithm, for K=3 factors.
% 
% Last updated: Sept. 18, 2011
% 

% global SCM;
% tol = 1e-8;


% initial estimates of X, Y, Z matrices
X_mat = eye(p);
Y_mat = eye(f);
Z_mat = eye(g);

% Likelihood = zeros(niter,1);
% Frob_error_ML_vec = zeros(niter,1);
% Frob_error_ML_cov_vec = zeros(niter,1);

% Theta_est = kron(X_mat,Y_mat); % initial guess

% temp1 = trace(Theta_est*SCM) - log(det(Theta_est));
% temp2 = norm(CovMat_inv-Theta_est,'fro')/norm(CovMat_inv,'fro');

% % compute ML solution (as in Werner's paper)
% Cov_inv_current = kron(X_mat,Y_mat);
% abs_diff = 10000;

for ni=1:niter,
%     Theta_prev = Theta_est;
    X_prev = X_mat;
    Y_prev = Y_mat;
    Z_prev = Z_mat;
        
    % compute temporary matrix T_X
    T_X = zeros(f,f);
    
    for i1=1:p,
        for i2=1:p,
            for j1=1:g,
                for j2=1:g,
                    T_X = T_X + X_mat(i2,i1)*Z_mat(j1,j2)*SCM(((i1-1)*f*g+1:g:i1*f*g) + j1-1,((i2-1)*f*g+1:g:i2*f*g) + j2-1);
                end
            end
        end
    end
    T_X = T_X/(p*g);
    B_hat = T_X;
    B_hat_inv = inv(T_X);    
    Y_mat = B_hat_inv;
    
    % compute T_Y
    T_Y = zeros(p,p);
    yrmert = kron(Y_mat,Z_mat);
    for j1=1:f*g,
        for j2=1:f*g,
        	T_Y = T_Y + yrmert(j2,j1)*SCM(j1:f*g:end,j2:f*g:end);
        end
    end
    T_Y = T_Y/(f*g);
    A_hat = T_Y;
    A_hat_inv = inv(T_Y);
    X_mat = A_hat_inv;
    
    % compute T_Z
    T_Z = zeros(g,g);
    zrmert = kron(X_mat,Y_mat);
    for i1=1:p*f,
        for i2=1:p*f,
        	T_Z = T_Z + zrmert(i2,i1)*SCM((i1-1)*g+1:i1*g,(i2-1)*g+1:i2*g);
        end
    end
    T_Z = T_Z/(p*f);
    C_hat = T_Z;
    C_hat_inv = inv(T_Z);
    Z_mat = C_hat_inv;
    
%     Theta_est = kron(X_mat, Y_mat);
    
%     Frob_diff = norm(Theta_prev-Theta_est,'fro');
  %  Frob_diff = sqrt(computeFrob(X_mat,Y_mat,X_prev,Y_prev));
    mt = kron(X_mat,kron(Y_mat,Z_mat));
    mtprv = kron(X_prev,kron(Y_prev,Z_prev));
    Frob_diff = norm(mt-mtprv,'fro');
    Frob_diff = Frob_diff./(norm(X_prev,'fro')*norm(Y_prev,'fro')*norm(Z_prev,'fro'));
    
%     Frob_diff, pause
    
    if Frob_diff<tol && ni > 2,
        disp(['Exited ML solution at FF iteration = ' num2str(ni)]);
%         pause
        break;
    end

%     % store objective value at end of each iteration and compute error in
%     % overall covariance matrix (to plot later)
%     Theta_est = kron(X_mat,Y_mat);

%     lik_est = trace(Theta_est*SCM) - log(det(Theta_est));
%     Likelihood(ni) = lik_est;
%     Frob_error_ML_vec(ni) = norm(CovMat_inv-Theta_est,'fro')^2/norm(CovMat_inv,'fro')^2; % inverse
%     Frob_error_ML_cov_vec(ni) = norm(CovMat-kron(A0,B_hat),'fro')^2/norm(CovMat,'fro')^2; % forward
    
end

% lik_optimal = trace(CovMat_inv*SCM) - log(det(CovMat_inv));
% Likelihood = [temp1; Likelihood] - lik_optimal*ones(niter+1,1);
% Frob_error_ML_vec = [temp2; Frob_error_ML_vec];

% Theta_0 = kron(X0,Y0);
% Sigma_0 = kron(A0,B0);
% Frob_error_ML_inv_final = norm(Theta_0-Theta_est,'fro')^2/norm(Theta_0,'fro')^2; % inverse
% Frob_error_ML_cov_final = norm(Sigma_0-kron(A_hat,B_hat),'fro')^2/norm(Sigma_0,'fro')^2; % forward

% whos
    
%temp1 = (norm(X0,'fro')*norm(Y0,'fro'))^2;
%temp2 = (norm(A0,'fro')*norm(B0,'fro'))^2;
%Frob_error_ML_inv_final = computeFrob(X_mat,Y_mat,X0,Y0)/temp1; % inverse
%Frob_error_ML_cov_final = computeFrob(A_hat,B_hat,A0,B0)/temp2; % forward


% Frob_error_ML_inv_final
% Frob_error_ML_cov_final
% pause



% % store output
% % ni_ML = ni;
% A_hat_ML = inv(X_mat);
% B_hat_ML = inv(Y_mat);
% CovMat_est_ML = kron(A_hat_ML, B_hat_ML);
% CovMat_est_inv_ML = Cov_inv_current;



A_mat = A_hat;
B_mat = B_hat;
C_mat = C_hat;


end

