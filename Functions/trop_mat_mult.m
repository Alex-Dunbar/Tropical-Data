function [out] = trop_mat_mult(A,x)
% Max-Plus matrix-vector multiplication
%
% Inputs: A - matrix of size m x n
%         x - vector of size n x 1
%
% Outputs: out - vector of size m x 1
%
% Example usage
%
% A = [0 1; 0 2]; x = [1;2];
% out = trop_mat_mult(A,x);
% returns out [3; 4]

%form proper rows
terms = A + x';

%take maximum
out = max(terms,[],2);

end
