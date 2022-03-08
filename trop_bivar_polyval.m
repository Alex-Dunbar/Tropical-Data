function [val] = trop_bivar_polyval(X,coeffs, d)
%Evaluate a bivariate Tropical polynomial
%
%Inputs: X      - N x 2 matrix of evaluation points
%        coeffs - size (d(1)+1)(d(2)+1) vector of coefficients. The order
%                 of the monomials is 1,x,2x,...,d(1)x,y,x+y, 2x + y, ..., 
%                                            d(1)x + d(2)y
%        d      - vector of two positive integers giving the degree in each
%                 variable
%
%Outputs: val   - vector of lenth N consisting of the evaluation of the
%                 polynomial at the points specified by X



inc = d(1)+1; %size of blocks
val = trop_polyval(X(:,1), [coeffs(2:inc); coeffs(1)]);
for t = 1:d(2) %each "power" of variable 2
    val = max(val, t*X(:,2) + trop_polyval(X(:,1), [coeffs((t*inc+2):(t+1)*inc); coeffs(t*inc + 1)]));
end

end