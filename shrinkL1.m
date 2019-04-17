function a_shr = shrinkL1( A,lambda )
%SHRINKL1 Summary of this function goes here
%   Detailed explanation goes here
a = A - diag(diag(A));
a_shr = sign(a).*max(0,abs(a) - lambda) + diag(diag(A));



end

