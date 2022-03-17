function [coeffs] = trop_mmae_fit(x,y,d)
% Computes coefficients of the optimal subsolution of 
% 
%             ||q(x) - y||_infty 
% 
% where q is a single variable tropical polynomial of degree d 
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

N = numel(x);
%X = [x*(1:d) zeros(size(x))];

%find optimal subfit
w = trop_poly_subfit(x,y,d);

%Compute shift for coefficients based on error in subfit

%Evaluate polynomial at data points without forming vandermonde
approx = trop_polyval(x,w);

%Shift by 1/2 infty norm error
mu = 0.5*norm(y-approx,'inf');

coeffs = w + mu;

end