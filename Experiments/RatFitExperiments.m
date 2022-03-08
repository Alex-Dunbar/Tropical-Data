%% Single Variable tropical rational function fitting tests
% Section 1: Fitting Tropical line Using alternating method
% Section 2: Fitting Tropical rational function using alternating method
% Section 3: Fitting sine curve using alternating method
% Section 4: Fitting sine curve using gradient approximation

%% Section 1: Piecewise linear function
x = linspace(-1,12,200); 
y = max(x-2,3);
x = x'; y = y'; 
y1 = y + 0.5*randn(size(y));

figure(1)
plot(x,y1,'k.')
hold on
for max_iter = [1 5 25]
    [num_coeffs, den_coeffs] = trop_rat_fit(x,y1,max_iter,2);
    fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);
    plot(x,fit);
end
title('Approximation of a Tropical Line')
legend('data','iter = 1', 'iter = 5', 'iter = 25','Location','northwest')
hold off

%% Section 2: Tropical Rational function
y = max(x+2,4) - max(x-3,1); y3 = y + 0.25*randn(size(y));

figure(2)
plot(x,y3,'k.')
hold on
for max_iter = [1 5 25]
    [num_coeffs, den_coeffs] = trop_rat_fit(x,y3,max_iter,15);
    fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);
    plot(x,fit);
end
title('Approximation of a Tropical Rational Function')
legend('data','iter = 1', 'iter = 5', 'iter = 25')
hold off


%% Section 3: Sine
y = sin(x); y2 = y + 0.1*randn(size(y)); 

figure(3) %plot first few iterations
plot(x,y2,'k.')
hold on
for max_iter = [1 5 25]
    [num_coeffs, den_coeffs] = trop_rat_fit(x,y2,max_iter,16);
    fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);
    plot(x,fit);
end
legend('data','iter = 1', 'iter = 5', 'iter = 25')
title('Degree 16 rational fit--Alternating Scheme')
hold off

figure(4)
plot(x,y2,'k.')
hold on

%Many iterations/degree 16
maxiter = 1250;
d = 16;

%compute fit
[num_coeffs, den_coeffs,err, L2_err,update_norm,l2_loss_grad] = trop_rat_fit(x,y2,maxiter,d);
fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);

plot(x,fit,'b');
title('Degree 16 rational fit--Alternating Scheme')
legend('data','degree 16 fit')
hold off

figure(5)
semilogy(1:maxiter,err,'k')
hold on
semilogy(1:maxiter,L2_err,'b')
semilogy(1:maxiter,update_norm,'r')
title('Error for Degree 16')
xlabel('Iteration')
ylabel('Error (log scale)')
legend('\infty norm','2 norm','norm of update')
hold off



%% Section 4: L2 norm minimization
x = linspace(-1,12,200); x = x';
y = sin(x); y2 = y + 0.1*randn(size(y));

maxiter = 5000;
d = 16; 

[num_coeffs,den_coeffs, err, grad] = trop_rat_2_fit(x,y,maxiter,d);
figure(6)
plot(x,y2,'k.')
hold on
fit = trop_polyval(x,num_coeffs) - trop_polyval(x,den_coeffs);
plot(x,fit,'b')
title('Degree 16 rational fit--Gradient approximation')
hold off

figure(7)
semilogy(1:numel(err), err)
hold on
semilogy(1:numel(grad),grad)
title('Error in fit via Gradient approximation')
legend('Square L2 error','Gradient Norm')
hold off








