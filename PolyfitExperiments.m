%Testing Tropical polynomial fitting. These examples are taken from 
%
% Maragos and Theodosis 2019 Tropical Geometry and Piecewise 
% Linear approximation of Curves and Surfaces on Weighted Lattices

x = linspace(-1,12,200); 
y = max(x-2,3);
x = x'; y = y'; 

%Gaussian Noise
y1 = y + 0.5*randn(size(y));

coeffs1 = trop_poly_subfit(x,y1,1);
fit1 = trop_polyval(x,coeffs1);

coeffs2 = trop_mmae_fit(x,y1,1);
fit2 = trop_polyval(x,coeffs2);

figure(1)
plot(x,y1,'k.')
hold on
plot(x,fit1,'r')
plot(x,fit2,'g')
legend('data','LS fit','MMAE fit','Location','northwest')
title('Gaussian Noise')
hold off

%Uniform Noise
y2 = y + -0.5 + rand(size(y));

coeffs1 = trop_poly_subfit(x,y2,1);
fit1 = trop_polyval(x,coeffs1);

coeffs2 = trop_mmae_fit(x,y2,1);
fit2 = trop_polyval(x,coeffs2);

figure(2)
plot(x,y2,'k.')
hold on
plot(x,fit1,'r')
plot(x,fit2,'g')
legend('data','LS fit','MMAE fit','Location','northwest')
title('Uniform Noise')
hold off


