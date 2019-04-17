function dy = MCP_derivative(A, a, rho)
%Function to evaluate derivative of (MCP penalty - L1) with parameter a > 0 and
%regularization parameter rho. Matricized with no offdiagonal penalty.
%
%Kristjan Greenewald, 8/4/18
%
%
A = A - diag(diag(A));
dy = (abs(A) < rho*a).*sign(A).*rho.*(-abs(A)/(rho*a));


%if abs(t) >= rho*a
%    dy = 0;
%else
%    dy = sign(t)*rho*(- abs(t)/(rho*a));
%end
return;