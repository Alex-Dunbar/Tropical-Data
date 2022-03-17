function [coeffs] = trop_poly_subfit(x,y,d)
% Computes coefficients of the optimal subsolution of 
% 
%      ||q(x) - y||_p  s.t.  q(x) <= y 
% 
% where q is a single variable tropical polynomial of degree d 
% via a min-plus matrix vector product. 
% 
% (See Maragos and Theodosis 2019 Tropical Geometry and Piecewise 
% Linear approximation of Curves and Surfaces on Weighted Lattices)
%
% Inputs: x - vector of independent data
%         y - vector of dependent data of size(x)
%         d - positive integer, maximum degree of tropical polynomial q.
%
% Outputs: coeffs - vector of length d+1 containing the coefficients of q.
%                   The minimizing polynomial is then 
%
%                   q(t) = max(t + coeffs(2), 2*t + coeffs(3), ...
%                                 d*t + coeffs(d+1), coeffs(1));
%
% Example Usage: 
%
% x = [0; 1; 2]; y = [0;0;1]; d= 2;
%
% coeffs = trop_poly_subfit(x,y,d);
% returns coeffs = [0;-1;-3]


%Compute min-plus matrix-vector product without forming "Vandermonde"
coeffs = zeros(d+1,1);
coeffs(1) = min(y);
for i = 2:(d+1)
    coeffs(i) = min(y-(i-1)*x);
end


end