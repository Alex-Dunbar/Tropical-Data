function [num_coeffs, den_coeffs,err,class_error,update_norm] = trop_bivar_rat_fit(X,y,max_iter,d)
% Uses an alternating method to find a tropical rational function
% approximation to given data
% 
% Inputs: X        - N x 2 matrix of independent data observations
%         y        - vector of dependent data of length n
%         max_iter - maximum number of iterations
%         d        - length 2 vector of positive integers giving the 
%                    maximum degree in each variable for both the numerator
%                    and denominator
%
% Outputs: num_coeffs  - vector of length (d(1)+1)(d(2)+1) containing the 
%                        coefficients of the numerator 
%          den_coeffs  - vector of length (d(1)+1)(d(2)+1) containing the 
%                        coefficients of the denominator
%          err         - infinity norm error at each iteration
%          class_err   - classification error in the case that y is binary
%                        0,1
%          update_norm - L2 norm of update
%
%Example Usage
%
% X = [0 0; 1 0; 2 4]; y = [0;0;1]; d= [2 2];
%
% [num_coeffs, den_coeffs] = trop_bivar_rat_fit(X,y,d);
% fit = trop_bivar_polyval(X,num_coeffs,d) - trop_bivar_polyval(X,den_coeffs,d)

%initialize
num_coeffs = -Inf*ones((d(1)+1)*(d(2)+1),1);
den_coeffs = -Inf*ones((d(1)+1)*(d(2)+1),1); den_coeffs(1) = mean(y);

err = zeros(max_iter,1);
class_error = zeros(max_iter,1);
update_norm = zeros(max_iter,1);

%Iterate
for k = 1:max_iter
    num_coeffs_old = num_coeffs; den_coeffs_old = den_coeffs;

    num_coeffs = bivar_trop_polyfit(X, trop_bivar_polyval(X,den_coeffs,d) + y,d);
    den_coeffs = bivar_trop_polyfit(X, trop_bivar_polyval(X,num_coeffs,d) - y,d);

    update_norm(k) = sqrt(norm(num_coeffs_old - num_coeffs)^2 + norm(den_coeffs_old - den_coeffs)^2 );
    
    %compute the fit and l^infty loss
    fit = trop_bivar_polyval(X,num_coeffs,d) - trop_bivar_polyval(X,den_coeffs,d);
    err(k) = norm( fit - y,"inf");
    
    %compute classification error in the case that $y$ is binary 
    class_error(k) = sum(abs((y - round(fit))))/size(X,1);

end

end