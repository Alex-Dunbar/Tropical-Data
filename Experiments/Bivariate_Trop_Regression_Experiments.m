%%Experiments with bivariate tropical functions
%
% Section 1: Approximating Tropical line with no noise
% Section 2: Approximating Tropical line with noise
% Section 3: Approximating a high degree tropical polynomial
% Section 4: Approximating a tropical rational function
% Section 5: Low degree approximation of peaks
% Section 6: High degree approximation of peaks
% Section 7: Binary Classification
% Section 8: Binary Classification on concentric circles w/ overlap

%% Section 1: Fitting to bivariate tropical Line
[X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [1,1]; %This fit should be exact, and is

g = @(X,Y) max(X+1,Y + 3); f = @(X,Y) max(g(X,Y),2);
Z = f(X,Y);y = reshape(Z,[size(Z,1)^2,1]);


[coeffs] = bivar_trop_polyfit(data,reshape(Z,[size(Z,1)^2,1]),d);

fit = trop_bivar_polyval(data,coeffs,d);

figure(1)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (5,5) Approximation of max(x+1,y+3,2), No noise')
hold off


%% Section 2: Fitting to bivariate tropical line with noise

[X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [5,5]; 

g = @(X,Y) max(X+1,Y + 3); f = @(X,Y) max(g(X,Y),2);
Z = f(X,Y); Z = Z + 0.1*randn(size(Z));
y = reshape(Z,[size(Z,1)^2,1]);


%Find and evaluate tropical rational approximation
[coeffs] = bivar_trop_polyfit(data,reshape(Z,[size(Z,1)^2,1]),d);
fit = trop_bivar_polyval(data,coeffs,d);


figure(2)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (5,5) Approximation of max(x+1,y+3,2), Noise')
hold off

%% Section 3: High degree polynomial 

[X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [5,5]; 

g = @(X,Y) max(5*X+ Y + 1,2*X + 3*Y + 3); f = @(X,Y) max(g(X,Y),2);
Z = f(X,Y); Z = Z + 0.1*randn(size(Z));
y = reshape(Z,[size(Z,1)^2,1]);

%Find and evaluate tropical rational approximation
[coeffs] = bivar_trop_polyfit(data,reshape(Z,[size(Z,1)^2,1]),d);
fit = trop_bivar_polyval(data,coeffs,d);


figure(3)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (5,5) Approximation of max(5x+ y+1,2x + 3y+3,2), noise')
hold off

%% Section 4: Rational function

[X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [5,5]; 

g = @(X,Y) max(5*X+ Y + 1,2*X + 3*Y + 3); 
h = @(X,Y) max(3*X + Y +2, X + 5*Y - 1);
f = @(X,Y) max(g(X,Y),2) - max(h(X,Y),1);
Z = f(X,Y); Z = Z + 0.1*randn(size(Z));
y = reshape(Z,[size(Z,1)^2,1]);

max_iter = 100;

%Find and Evaluate Tropical Rational Approximation
[num_coeffs, den_coeffs,err,L2_err,update_norm] = trop_bivar_rat_fit(data,y,max_iter,d);
fit = trop_bivar_polyval(data,num_coeffs,d) - trop_bivar_polyval(data,den_coeffs,d);

figure(4)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (5,5) Approximation of max(5x+ y+1,2x + 3y+3,2), noise')
hold off

figure(5)
semilogy(1:max_iter,err)
hold on 
semilogy(1:max_iter,update_norm)
legend('error','update norm')
title('Error in approximation of noisy tropical rational function')
hold off

%% Section 5: Fitting Rational Function to Peaks data--degree (10,10);

[X,Y,Z] = peaks;
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [10,10];
y = reshape(Z,[size(Z,1)^2,1]);

max_iter = 100;

[num_coeffs, den_coeffs,err,~,update_norm] = trop_bivar_rat_fit(data,y,max_iter,d);
fit = trop_bivar_polyval(data,num_coeffs,d) - trop_bivar_polyval(data,den_coeffs,d);


figure(6)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (10,10) Approximation of Peaks')
hold off

figure(7)
semilogy(1:max_iter,err)
hold on 
semilogy(1:max_iter,update_norm)
legend('error','update norm')
title('Error in degree (10,10) approximation of Peaks')
hold off

%% Section 6: Fitting Rational Function to Peaks data--degree (35,35);

[X,Y,Z] = peaks;
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [35,35];
y = reshape(Z,[size(Z,1)^2,1]);

max_iter = 120;

[num_coeffs, den_coeffs,err,~,update_norm] = trop_bivar_rat_fit(data,y,max_iter,d);
fit = trop_bivar_polyval(data,num_coeffs,d) - trop_bivar_polyval(data,den_coeffs,d);


figure(8)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (35,35) Approximation of Peaks')
hold off

figure(9)
semilogy(1:max_iter,err)
hold on 
semilogy(1:max_iter,update_norm)
legend('error','update norm')
title('Error in degree (35,35) approximation of Peaks')
hold off

%% Section 7: Bivariate Classification
X1 = ones(100,2) + 0.5*randn(100,2);
X2 = 0.5*randn(100,2);
data = [X1; X2];
y = [ones(100,1); zeros(100,1)];

d = [5,5];
max_iter = 100;

% Find coefficients
[num_coeffs, den_coeffs,~,class_error,update_norm] = trop_bivar_rat_fit(data,y,max_iter,d);

%Find training Error
fit = trop_bivar_polyval(data,num_coeffs,d) - trop_bivar_polyval(data,den_coeffs,d);
predictions = round(fit);
train_error = sum(abs(y - predictions))/200


%Evaluate on new data
[X_plot,Y_plot] = meshgrid([-2:0.1:2]',[-2:0.1:2]');
plot_data = [reshape(X_plot,[size(X_plot,1)^2,1]) reshape(Y_plot,[size(Y_plot,1)^2,1])];
fit = trop_bivar_polyval(plot_data,num_coeffs,d) - trop_bivar_polyval(plot_data,den_coeffs,d);


figure(10)
surf(X_plot,Y_plot,reshape(fit,[size(X_plot,1) size(X_plot,2)]));
hold on
plot3(data(1:100,1),data(1:100,2),y(1:100),'b.');
plot3(data(101:200,1),data(101:200,2),y(101:200),'r.');
legend('Approximation','Class 1','Class 2')
title('Degree (5,5) Approximation')
hold off

figure(11)
semilogy(1:max_iter,class_error)
hold on 
semilogy(1:max_iter,update_norm)
legend('classification error','update norm')
title('Error in degree (5,5) approximation')
hold off


%% Section 8: Bivariate Classification: Concentric Circles

N = 100;
thetas = 2*pi*rand(2*N,1);
r1 = 1.1*rand(N,1); r2 = 1 + 0.5*rand(N,1);

X1 = [r1.*cos(thetas(1:N))   r1.*sin(thetas(1:N))]; 
X2 = [r2.*cos(thetas((N+1):2*N)) r2.*sin(thetas((N+1):2*N))];


data = [X1; X2];
y = [ones(N,1); zeros(N,1)];

d = [5,5]; 
max_iter = 100;

% Find coefficients
[num_coeffs, den_coeffs,~,class_error,update_norm] = trop_bivar_rat_fit(data,y,max_iter,d);

%Find training Error
fit = trop_bivar_polyval(data,num_coeffs,d) - trop_bivar_polyval(data,den_coeffs,d);
predictions = round(fit);
train_error = sum(abs(y - predictions))/(2*N)


%Evaluate on new data
[X_plot,Y_plot] = meshgrid([-2:0.1:2]',[-2:0.1:2]');
plot_data = [reshape(X_plot,[size(X_plot,1)^2,1]) reshape(Y_plot,[size(Y_plot,1)^2,1])];
fit = trop_bivar_polyval(plot_data,num_coeffs,d) - trop_bivar_polyval(plot_data,den_coeffs,d);


figure(12)
surf(X_plot,Y_plot,reshape(fit,[size(X_plot,1) size(X_plot,2)]));
hold on
plot3(data(1:N,1),data(1:N,2),y(1:N),'b.');
plot3(data((N+1):2*N,1),data((N+1):2*N,2),y((N+1):(2*N)),'r.');
legend('Approximation','Class 1','Class 2')
title('Degree (5,5) Approximation')
hold off

figure(13)
semilogy(1:max_iter,class_error)
hold on 
semilogy(1:max_iter,update_norm)
legend('classification error','update norm')
title('Error in degree (5,5) approximation')
hold off

figure(14)
plot(X1(:,1),X1(:,2),'b.')
hold on 
plot(X2(:,1),X2(:,2),'r.')
title('Training Data')
hold off





