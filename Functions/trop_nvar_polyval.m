function [val] = trop_nvar_polyval(X,coeffs, d)
%Evaluate a bivariate Tropical polynomial
%
%Inputs: X      - N x n matrix of evaluation points
%        coeffs - size (d(1)+1)(d(2)+1)...(d(n)+1) vector of coefficients. The order
%                 of the monomials is 1,x,2x,...,d(1)x,y,x+y, 2x + y, ..., 
%                                            d(1)x + d(2)y, ...
%        d      - vector of n positive integers giving the degree in each
%                 variable
%
%Outputs: val   - vector of lenth N consisting of the evaluation of the
%                 polynomial at the points specified by X
%
% Example Usage: 
%
% X = [0 0; 1 0; 2 4]; y = [0;0;1]; d= [2 2];
%
% coeffs = trop_nvar_polyfit(X,y,d);
% fit = trop_nvar_polyval(X,coeffs,d)


n = numel(d); N = size(X,1);
if n == 1
    deg = numel(coeffs) - 1;
    val = coeffs(1)*ones(N,1);
    for t = 1:d
        val = max(val,t*X(:,1) + coeffs(t+1));
    end
else

inc = prod(d(1:n-1)+1); %size of blocks

val = trop_nvar_polyval(X(:,1:n-1), coeffs(1:inc),d(1:n-1));
for t = 1:d(n) %each "power" of variable n
    val = max(val, t*X(:,n) + trop_nvar_polyval(X(:,1:n-1), coeffs((t*inc+1):(t+1)*inc),d(1:n-1)));
end

end