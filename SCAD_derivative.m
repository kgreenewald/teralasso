function dy = SCAD_derivative(A, a, rho)
%Function to evaluate derivative of (SCAD penalty - L1) with parameter a > 2 and
%regularization parameter rho. Matricized with no offdiagonal penalty.
%
%Kristjan Greenewald, 8/4/18
%
%
A = A - diag(diag(A));
trm1 = (abs(A) > rho & abs(A) <= a*rho).*(-2*A./(2*(a-1)) + rho.*sign(A)/(a-1));
trm2 = (abs(A) > a*rho).*(-rho.*sign(A));
dy = trm1 + trm2;

%if abs(t) <= rho
%    dy = 0;
%elseif abs(t) <= a*rho
%    dy = -2*t/(2*(a-1)) + rho*sign(t)/(a-1);
%else
%    dy = - rho*sign(t);
%end
return;