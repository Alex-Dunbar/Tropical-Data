function [coeffs,eval] = trop_nvar_polyfit(X,y,d)
% Computes coefficients of the optimal solution of 
% 
%               ||q(x) - y||_infty
% 
% where q is a multivariate variable tropical polynomial of degree d 
% via min-plus matrix vector product. 
% 
% (See Maragos and Theodosis 2019 Tropical Geometry and Piecewise 
% Linear approximation of Curves and Surfaces on Weighted Lattices)
%
% Inputs: X - N x n matrix of independent data observations
%         y - vector of dependent data of length N
%         d - length 2 vector of positive integers giving the maximum 
%             degree in each variable for the tropical polynomial q.
%
% Outputs: coeffs - vector of length (d(1)+1)(d(2)+1)...(d(n)+1) containing the 
%                   coefficients of q.
%
%                   The minimizing polynomial is then 
%
%                   q(t,s,...) = max(coeffs(1), t + coeffs(2), ... , 
%                                t+s+ coeffs(d(1)+3), ..., 
%                                d(1)t + d(2)s + coeffs((d(1)+1)(d(2)+1), ...,);
%
% Example Usage: 
%
% X = [0 0; 1 0; 2 4]; y = [0;0;1]; d= [2 2];
%
% coeffs = trop_nvar_polyfit(X,y,d);
% fit = trop_nvar_polyval(X,coeffs,d)

%make sure d is column vector
d = reshape(d,[numel(d),1]);

row = 0; k = prod(d+1);
p = zeros(k,1);
n = numel(d); last = zeros(n,1);
eval = zeros(size(X,1),1);

%Compute min-plus matrix-vector product without forming "Vandermonde"

while any(last < d)
    u = X*last;
    for j = 1:n
        if last(j)<d(j)
            row = row + 1;
            %w(row) = min(y - X*last);
            p(row) = min(y - u);
            if row == 1
                eval = p(row)*ones(size(eval));
            else
                eval = max(eval, p(row) + u); %update evaluation of polynomial
            end
            last(j) = last(j)+ 1;
            u = u + X(j); %updated (X*last)
            if j > 1
                last(1:j-1) = 0;
                break
            end
            break
        end
    end
end

row = row + 1;
p(row) = min(y - X*d);
eval = max(eval,p(row)+X*d);

%shift to minimize infty norm
mu = 0.5*norm(y-eval,"inf");
coeffs = p + mu;
eval = eval + mu;

end
