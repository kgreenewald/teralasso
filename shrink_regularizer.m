function a_shr = shrink_regularizer(A, lambda,stepSz, type,a)
%type = L1, SCAD, MCP

switch type
    case 'L1'
        a_shr = shrinkL1(A, lambda*stepSz);
      
    case 'SCAD'
        a_shr = shrinkL1(A - stepSz*SCAD_derivative(A,a,lambda), lambda*stepSz);
    case 'MCP'
        a_shr = shrinkL1(A - stepSz*MCP_derivative(A,a,lambda), lambda*stepSz);
end