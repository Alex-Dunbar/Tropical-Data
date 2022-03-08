function [num_grad,den_grad] = trop_rat_2_loss_grad(num_coeffs,den_coeffs,x,y)
% 
% Approximates the gradient of 
% 
%            L = sum_j (p(x_j) - q(x_j) - y_j)^2, 
% 
% where p,q are univariate tropical polynomials specified by num_coeffs 
% and den_coeffs, respectively.
%
% Inputs: num_coeffs - vector of coeffcients for the numerator 
%                      (i.e. p(x_j) = trop_polyval(x_j,num_coeffs))
%         den_coeffs - vector of coefficients for the denominator
%                      (i.e. q(x_j) = trop_polyval(x_j,den_coeffs))
%         x          - vector of independent variable data
%         y          - vector of dependent variable data
%
% Outputs: num_grad  - vector of size(num_coeffs) giving the gradient of
%                      L with respect to numerator coefficients
%          den_grad  - vector of size(den_coeffs) giving the gradient of L
%                      with respect to denominator coefficients


N = numel(y); d = numel(num_coeffs) -1;

num_vec = trop_polyval(x,num_coeffs);
den_vec = trop_polyval(x,den_coeffs);
magnitude = 2*(num_vec- den_vec - y);

%Find where maximum is achieved in numerator and denominator-probably a
%better way to do this than brute force search. 
num_indicator = zeros(d+1,N);
for j = 1:N
    for ell = 1:d
        num_indicator(ell,j) = (ell*x(j) + num_coeffs(ell)) == num_vec(j);
    end
    num_indicator(d+1,j) = num_coeffs(d+1)==num_vec(j);
    %scale
    num_indicator(:,j) = num_indicator(:,j)/sum(num_indicator(:,j));
end

num_indicator = sparse(num_indicator);

den_indicator = zeros(d+1,N);
for j = 1:N
    for ell = 1:d
        den_indicator(ell,j) = (ell*x(j) + den_coeffs(ell)) == den_vec(j);
    end
    den_indicator(d+1,j) = den_coeffs(d+1)==den_vec(j);
    %scale
    den_indicator(:,j) = den_indicator(:,j)/sum(den_indicator(:,j));
end
den_indicator = sparse(den_indicator);

num_grad = num_indicator*magnitude;
den_grad = (-1)*den_indicator*magnitude;

end
