function [num_coeffs, den_coeffs,err,L2_err,update_norm,l2_loss_grad] = trop_rat_fit(x,y,max_iter,d)
% Fits a tropical rational function to the data x,y using an alternating
% method to find the numerator and denominator coefficients
%
% Inputs: x        - vector of independent data
%         y        - vector of dependent data of size(x)
%         max_iter - integer, maximum number of iterations
%         d        - integer, maximum degree for numerator and denominator
%
% Outputs: num_coeffs   - vector of length d+1 with coefficients for the 
%                         numerator of the rational fit 
%          den_coeffs   - vector of length d+1 with coefficients for the
%                         denominator of the rational fit
%          err          - vector with infinity norm error at each iteration
%          L2_err       - vector with L2 error at each iteration
%          update_norm  - vector with L2 norm of update at each iteration
%          l2_loss_grad - approximate gradient of the L2 loss at each
%                        iteration
%
%Example usage:
%
% x = linspace(1,100); x = x'; y = sin(x);
% [num_coeffs, den_coeffs] = trop_rat_fit(x,y,100,16)
% 
% fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs)
%

%initialize coefficients
num_coeffs = -Inf*ones(d+1,1);den_coeffs = [-Inf*ones(d,1);-mean(y)];

%initialize error outputs
err = zeros(max_iter,1);
L2_err = zeros(max_iter,1);
update_norm = zeros(max_iter,1);
l2_loss_grad = zeros(max_iter,1);

%Iterate fitting numerator and denominator
for k = 1:max_iter
    num_coeffs_old = num_coeffs; den_coeffs_old = den_coeffs;
    num_coeffs = trop_mmae_fit(x, trop_polyval(x,den_coeffs) + y,d);
    den_coeffs = trop_mmae_fit(x, trop_polyval(x,num_coeffs) - y,d);

    update_norm(k) = sqrt(norm(num_coeffs_old - num_coeffs)^2 + norm(den_coeffs_old - den_coeffs)^2 );
    
    %compute the fit, l^infty loss, and l^2 loss
    fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);
    err(k) = norm( fit - y,"inf");
    L2_err(k) = norm(fit - y)/norm(y);

    [num_grad,den_grad] = trop_rat_2_loss_grad(num_coeffs,den_coeffs,x,y);
    l2_loss_grad(k) = sqrt(norm(num_grad)^2 + norm(den_grad)^2);
end

end
