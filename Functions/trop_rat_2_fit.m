function [num_coeffs,den_coeffs, err, grad] = trop_rat_2_fit(x,y,maxiter,d)
% 
% Fits a tropical rational function to the data x,y using an approximation
% to the gradient of 
%
%            L = sum_j (p(x_j) - q(x_j) - y_j)^2, 
% 
%
% Inputs: x        - vector of independent data
%         y        - vector of dependent data of size(x)
%         maxiter  - integer, maximum number of iterations
%         d        - integer, maximum degree for numerator and denominator
%
% Outputs: num_coeffs   - vector of length d+1 with coefficients for the 
%                         numerator of the rational fit 
%          den_coeffs   - vector of length d+1 with coefficients for the
%                         denominator of the rational fit
%          err          - vector with squared L2 norm of loss at each
%                         iteration
%          grad         - vector with L2 norm of gradient approximation at
%                         each iteration
%
%Example usage:
%
% x = linspace(1,100); x = x'; y = sin(x);
% [num_coeffs, den_coeffs] = trop_rat_2_fit(x,y,100,16)
% 
% fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs)
%

%initialize -Initialization seems to be very important here
[num_coeffs,den_coeffs] = trop_rat_fit(x,y,25,d);

err = zeros(maxiter,1);
grad = zeros(maxiter,1);
grad(1) = 100;

iter = 1;
while grad(max(iter-1,1)) >= 10^(-5) && iter < maxiter
    [num_grad,den_grad] = trop_rat_2_loss_grad(num_coeffs,den_coeffs,x,y);

    %iterate
    num_coeffs = num_coeffs - 0.001*num_grad; den_coeffs = den_coeffs - 0.001*den_grad;

    err(iter) = norm(trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs) - y)^2;
    grad(iter) = sqrt(norm(num_grad)^2 + norm(den_grad)^2);
    iter = iter + 1;
end

err = err(err > 0);
grad = grad(grad > 0);

end
