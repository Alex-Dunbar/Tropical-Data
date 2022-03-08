function [coeffs] = bivar_trop_polyfit(X,y,d)
% Computes coefficients of the optimal subsolution of 
% 
%               ||q(x) - y||_infty
% 
% where q is a bivariate variable tropical polynomial of degree d 
% via a min-plus matrix vector product. 
% 
% (See Maragos and Theodosis 2019 Tropical Geometry and Piecewise 
% Linear approximation of Curves and Surfaces on Weighted Lattices)
%
% Inputs: X - N x 2 matrix of independent data observations
%         y - vector of dependent data of length n
%         d - length 2 vector of positive integers giving the maximum 
%             degree in each variable for the tropical polynomial q.
%
% Outputs: coeffs - vector of length (d(1)+1)(d(2)+1) containing the 
%                   coefficients of q.
%
%                   The minimizing polynomial is then 
%
%                   q(t,s) = max(coeffs(1), t + coeffs(2), ... , 
%                           t+s+ coeffs(d(1)+3), ..., 
%                           d(1)t + d(2)s + coeffs((d(1)+1)(d(2)+1));
%
% Example Usage: 
%
% X = [0 0; 1 0; 2 4]; y = [0;0;1]; d= [2 2];
%
% coeffs = bivar_torp_polyfit(X,y,d);
% fit = trop_bivar_polyval(X,coeffs,d)



row = 0; k = (d(1)+1)*(d(2)+1);
w = zeros(k,1);


%Compute min-plus matrix-vector product without forming "Vandermonde"
for t = 0:d(2)
    for s = 0:d(1)
        row = row+1;
        w(row) = min(y-s*X(:,1)-t*X(:,2));
    end
end

%shift to minimize infty norm
mu = 0.5*norm(trop_bivar_polyval(X,w,d) - y,'inf');
coeffs = w + mu;

end
