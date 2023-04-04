function [num_coeffs, den_coeffs, prob_out] = trop_nvar_rat_fit(X,y,max_iter,d,tol,options,test_data,test_y)
% Uses an alternating method to find a tropical rational function
% approximation to given data
% 
% Inputs: X        - N x n matrix of independent data observations
%         y        - vector of dependent data of length N
%         max_iter - maximum number of iterations
%         d        - length n vector of positive integers giving the 
%                    maximum degree in each variable for both the numerator
%                    and denominator
%         tol      - tolerance for stopping criterion based on update norm
%         options  - struct describing needed output computations:
%                     err - infinity norm error at each iteration
%                     update - update norm at each iteration
%                     test- error on test set at each iteration
%                     class - classification error at each iteration
%                     L2 - L2 error at each iteration
%         tol      - stopping criterion tolerance for update norm
%         test_data- n column matrix of test data independent values
%         test_y   - vector of test data response values.
%
% Outputs: num_coeffs  - vector of length (d(1)+1)(d(2)+1)...(d(n)+1) containing the 
%                        coefficients of the numerator 
%          den_coeffs  - vector of length (d(1)+1)(d(2)+1)...(d(n)+1) containing the 
%                        coefficients of the denominator
%          
%          prob_out  - struct containing the following fields:
%          err         - infinity norm error at each iteration
%          update      - infinity norm of update
%          test        - error on test set at each iteration
%          class       - classification error in the case that y is binary
%                        0,1
%          L2          - L2 error at each iteration
%          iter        - number of iterations performed
%
%Example Usage
%
% [X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
% data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
% d = [5,5]; 
% 
% max_iter = 100;
% options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
% tol = 10^(-12);
%
% [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);
% fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);

n = size(X,2);N = size(X,1);

%initialize
num_coeffs = -Inf*ones(prod(d + 1),1);
den_coeffs = -Inf*ones(prod(d + 1),1); den_coeffs(1) = -mean(y);

prob_out = struct();
if options.err == 1
    prob_out.err = zeros(max_iter,1);
end
if options.update == 1
    prob_out.update = zeros(max_iter,1); 
end
if options.test == 1
    prob_out.test = zeros(max_iter,1);
end
if options.class == 1
    prob_out.class = zeros(max_iter,1);
end
if options.L2 == 1
    prob_out.L2 = zeros(max_iter,1);
end

num_val = trop_nvar_polyval(X,num_coeffs,d);
den_val = trop_nvar_polyval(X,den_coeffs,d);

%Iterate
for k = 1:max_iter
    num_coeffs_old = num_coeffs; den_coeffs_old = den_coeffs;
    
    num_coeffs = trop_nvar_polyfit(X, den_val + y,d);
    num_val = trop_nvar_polyval(X,num_coeffs,d);

    den_coeffs = trop_nvar_polyfit(X, num_val - y,d);
    den_val = trop_nvar_polyval(X,den_coeffs,d);
    
    update_norm = max(norm(num_coeffs_old - num_coeffs, "inf"), norm(den_coeffs_old-den_coeffs,"inf"));
    
    %store update norm if needed
    if options.update == 1
        prob_out.update(k) = update_norm;
    end

    %compute the fit and l^infty loss
    if options.err == 1
        fit = num_val - den_val;
        prob_out.err(k) = norm( fit - y,"inf");
    end
    
    %L2 loss
    if options.L2 == 1
        prob_out.L2(k) = norm( fit - y);
    end

    %compute and store test error if needed
    if options.test == 1
        test_fit = trop_nvar_polyval(test_data,num_coeffs,d) - trop_nvar_polyval(test_data,den_coeffs,d);
        prob_out.test(k) = norm(test_fit - test_y, "inf");
    end
    
    %compute classification error in the case that y is binary 
    if options.class == 1
        fit = num_val - den_val;
        prob_out.class(k) = sum(abs((y - round(fit))))/size(X,1);
    end
    
    %stop iterations at tolerance
    if update_norm < tol
        prob_out.iterations = k;
        return
    end
end

prob_out.iterations = max_iter;

end