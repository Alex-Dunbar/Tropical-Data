%%Experiments with multivariate tropical functions
%
% Section 1: Approximating Bivariate Tropical line with no noise
% Section 2: Approximating Bivariate Tropical line with noise
% Section 3: Approximating a high degree bivariate tropical polynomial
% Section 4: Approximating a bivariate tropical rational function
% Section 5: Low degree approximation of peaks
% Section 6: High degree approximation of peaks
% Section 7: Binary Classification
% Section 8: Binary Classification on concentric circles w/ overlap
% Section 9: Recovery of Univariate and Bivariate Tropical Rational Functions
% Section 10: Recovery of 6 variable Tropical Rational Functions
% Section 11: Recovery of 6 variable Tropoical Rational Functions (with noise)
% Section 12: More Bivariate Examples

%% Section 1: Fitting to bivariate tropical Line
[X,Y] = meshgrid([-3:0.125:3]',[-3:0.125:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [1,1]; %This fit should be exact, and is

g = @(X,Y) max(X+1,Y + 3); f = @(X,Y) max(g(X,Y),2);
Z = f(X,Y);y = reshape(Z,[size(Z,1)^2,1]);

[coeffs] = trop_nvar_polyfit(data,y,d);
fit = trop_nvar_polyval(data,coeffs,d);

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
[coeffs] = trop_nvar_polyfit(data,y,d);
fit = trop_nvar_polyval(data,coeffs,d);


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
[coeffs] = trop_nvar_polyfit(data,y,d);
fit = trop_nvar_polyval(data,coeffs,d);

figure(3)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
legend('Approximation','data')
title('Degree (5,5) Approximation of max(5x+ y+1,2x + 3y+3,2), noise')
hold off

%% Section 4: Rational function

[X,Y] = meshgrid([-3:0.1:3]',[-3:0.1:3]');
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [5,5]; 

g = @(X,Y) max(5*X+ Y + 1,2*X + 3*Y + 3); 
h = @(X,Y) max(3*X + Y +2, X + 5*Y - 1);
f = @(X,Y) max(g(X,Y),2) - max(h(X,Y),1);
Z = f(X,Y); %Z = Z  + 0.1*randn(size(Z));
y = reshape(Z,[size(Z,1)^2,1]);

max_iter = 100;
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);

%Find and Evaluate Tropical Rational Approximation
[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);
fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);

figure(4)
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
hold on
plot3(data(:,1),data(:,2),y,'k.');
%legend('Approximation','data')
%title('Degree (5,5) Approximation of max(5x+ y+1,2x + 3y+3,2), noise')
%matlab2tikz('RationalFunctionExample.tex','showInfo',false)
hold off

figure(5)
semilogy(1:max_iter,prob_out.err)
hold on 
semilogy(1:max_iter,prob_out.update)
legend('error','update norm')
title('Error in approximation of noisy tropical rational function')
hold off

%% Section 5: Fitting Rational Function to Peaks data--degree (10,10);

[X,Y,Z] = peaks;
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [10,10];
y = reshape(Z,[size(Z,1)^2,1]);
%y = y+0.25*randn(size(y)); %uncomment to add noise to data

max_iter = 250;
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);

[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);
fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);

%uncomment to get parameters for expression as NN
%writematrix([0;num_coeffs],"peaks_num.csv")
%writematrix([0;den_coeffs],"peaks_den.csv")

%Exponent_mat = [0 0];
%for j = 0:10
%    for i = 0:10
%        Exponent_mat = [Exponent_mat; i j];
%    end
%end

%writematrix(Exponent_mat,"peaks_ex.csv")
        

fig6 = figure(6);
plot3(data(:,1),data(:,2),y,'k.');
hold on
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
legend('Data','Approximation')
title('Degree 10 Tropical Rational Fit')


set(fig6,'Units','Inches');
pos = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig6,'peaks_1010_alternate_fit','-dpdf','-r0')
hold off
matlab2tikz('tikz/peaks_1010_alternate_fit.tex','showInfo',false)

fig7 = figure(7);
semilogy(1:prob_out.iterations,prob_out.err(1:prob_out.iterations),'k')
hold on 
semilogy(1:prob_out.iterations,prob_out.update(1:prob_out.iterations),'r')
xlabel('Iteration')
ylabel('Error (log scale)')
legend('\infty-norm error','update norm','Location','southwest')
title('Error for Degree 10 Fit to Peaks Data')
% set(fig7,'Units','Inches');
% pos = get(fig7,'Position');
% set(fig7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig7,'peaks_1010_alternate_convergence','-dpdf','-r0')

hold off
matlab2tikz('tikz/peaks_1010_alternate_convergence.tex','showInfo',false)

%% Section 6: Fitting Rational Function to Peaks data--degree (31,31);

[X,Y,Z] = peaks;
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [31,31];
y = reshape(Z,[size(Z,1)^2,1]);
%y = y+0.25*randn(size(y)); %Uncomment to add noise

max_iter = 500;
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);

[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);
fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);

%uncomment to get parameters for expression as NN
%writematrix([0;num_coeffs],"peaks_num.csv")
%writematrix([0;den_coeffs],"peaks_den.csv")

%Exponent_mat = [0 0];
%for j = 0:31
%    for i = 0:31
%        Exponent_mat = [Exponent_mat; i j];
%    end
%end

%writematrix(Exponent_mat,"peaks_ex.csv")
        
fig6 = figure(6);
plot3(data(:,1),data(:,2),y,'k.');
hold on
surf(X,Y,reshape(fit,[size(X,1) size(X,2)]));
legend('Data','Approximation')
title('Degree 31 Tropical Rational Fit')

set(fig6,'Units','Inches');
pos = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig6,'peaks_3131_alternate_fit','-dpdf','-r0')
hold off
matlab2tikz('tikz/peaks_3131_alternate_fit.tex','showInfo',false)

fig7 = figure(7);
semilogy(1:max_iter,prob_out.err,'k')
hold on 
semilogy(1:max_iter,prob_out.update,'r')
xlabel('Iteration')
ylabel('Error (log scale)')
legend('\infty-norm error','update norm','Location','southwest')
title('Error for Degree 31 Fit to Peaks Data')
set(fig7,'Units','Inches');
pos = get(fig7,'Position');
set(fig7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig7,'peaks_3131_alternate_convergence','-dpdf','-r0')
hold off
matlab2tikz('tikz/peaks_3131_alternate_convergence.tex','showInfo',false)

%% Section 7: Fitting Rational Function to Peaks data--degree (35,35) & trying scaling parameters;

[X,Y,Z] = peaks;
data = [reshape(X,[size(X,1)^2,1]) reshape(Y,[size(Y,1)^2,1])];
d = [35,35];
y = reshape(Z,[size(Z,1)^2,1]); 

max_iter = 5000;
options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
tol = 10^(-12);

kvec = 1:0.2:3;
err_k = zeros(numel(kvec),1);
iter = 0;
for k = kvec
    disp(k)
    iter = iter + 1;
    [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(k*data,y,max_iter,d,tol,options);
    fit = trop_nvar_polyval(k*data,num_coeffs,d) - trop_nvar_polyval(k*data,den_coeffs,d);
    
    err_k(iter) = norm(fit - y, "inf");
end

fig8 = figure(8);
plot(kvec,err_k)
hold on
title('Error vs Scaling Parameter')
ylabel('Error')
xlabel('c')
legend('Error')

% set(fig8,'Units','Inches');
% pos = get(fig8,'Position');
% set(fig8,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig8,'scaling parameter','-dpdf','-r0')
% hold off
% matlab2tikz('tikz/scaling parameter.tex','showInfo',false)

%% Section 8: Bivariate Classification
X1 = ones(100,2) + 0.5*randn(100,2);
X2 = 0.5*randn(100,2);
data = [X1; X2];
y = [ones(100,1); zeros(100,1)];

d = [5,5];
max_iter = 100;
options = struct('err',0,'update',1,'test',0,'class',1,'L2',0);
tol = 10^(-12);

% Find coefficients
[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);

%Find training Error
fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);
predictions = round(fit);
train_error = sum(abs(y - predictions))/200


%Evaluate on new data
[X_plot,Y_plot] = meshgrid([-2:0.1:2]',[-2:0.1:2]');
plot_data = [reshape(X_plot,[size(X_plot,1)^2,1]) reshape(Y_plot,[size(Y_plot,1)^2,1])];
fit = trop_nvar_polyval(plot_data,num_coeffs,d) - trop_nvar_polyval(plot_data,den_coeffs,d);


figure(10)
surf(X_plot,Y_plot,reshape(fit,[size(X_plot,1) size(X_plot,2)]));
hold on
plot3(data(1:100,1),data(1:100,2),y(1:100),'b.');
plot3(data(101:200,1),data(101:200,2),y(101:200),'r.');
legend('Approximation','Class 1','Class 2')
title('Degree (5,5) Approximation')
hold off

figure(11)
semilogy(1:max_iter,prob_out.class(1:prob_out.iterations))
hold on 
semilogy(1:max_iter,prob_out.update(1:prob_out.iterations))
legend('classification error','update norm')
title('Error in degree (5,5) approximation')
hold off


%% Section 9: Bivariate Classification: Concentric Circles

N = 100;
thetas = 2*pi*rand(2*N,1);
r1 = 1.1*rand(N,1); r2 = 1 + 0.5*rand(N,1);

X1 = [r1.*cos(thetas(1:N))   r1.*sin(thetas(1:N))]; 
X2 = [r2.*cos(thetas((N+1):2*N)) r2.*sin(thetas((N+1):2*N))];


data = [X1; X2];
y = [ones(N,1); zeros(N,1)];

d = [5,5]; 
max_iter = 100;
options = struct('err',0,'update',1,'test',0,'class',1,'L2',0);
tol = 10^(-12);


% Find coefficients
[num_coeffs, den_coeffs, prob_out] = trop_nvar_rat_fit(data,y,max_iter,d,tol,options);

%Find training Error
fit = trop_nvar_polyval(data,num_coeffs,d) - trop_nvar_polyval(data,den_coeffs,d);
predictions = round(fit);
train_error = sum(abs(y - predictions))/(2*N)


%Evaluate on new data
[X_plot,Y_plot] = meshgrid([-2:0.1:2]',[-2:0.1:2]');
plot_data = [reshape(X_plot,[size(X_plot,1)^2,1]) reshape(Y_plot,[size(Y_plot,1)^2,1])];
fit = trop_nvar_polyval(plot_data,num_coeffs,d) - trop_nvar_polyval(plot_data,den_coeffs,d);


figure(12)
surf(X_plot,Y_plot,reshape(fit,[size(X_plot,1) size(X_plot,2)]));
hold on
plot3(data(1:N,1),data(1:N,2),y(1:N),'b.');
plot3(data((N+1):2*N,1),data((N+1):2*N,2),y((N+1):(2*N)),'r.');
legend('Approximation','Class 1','Class 2')
title('Degree (5,5) Approximation')
hold off

figure(13)
semilogy(1:max_iter, prob_out.class(1:prob_out.iterations))
hold on 
semilogy(1:max_iter, prob_out.update(1:prob_out.iterations))
legend('classification error','update norm')
title('Error in degree (5,5) approximation')
hold off

figure(14)
plot(X1(:,1),X1(:,2),'b.')
hold on 
plot(X2(:,1),X2(:,2),'r.')
title('Training Data')
hold off

%% 6-variable problem
nvar = 6;
dmax = 3;

rng(1243093)%for reproducibility 6-var

f = @(x) x(1)*x(2)*x(3) + 2*x(4)*x(5)^2*sin(x(6)); %- exp(x(7)*x(8)*x(9)*x(10));
d = dmax*ones(nvar,1);

N = 10000; X = rand(N,nvar); 

train_data = zeros(N,1);
for index = 1:N
    train_data(index) = f(X(index,:));
end

test_X = rand(N,nvar); 

test_y = zeros(N,1);
for index = 1:N
   test_y(index) = f(test_X(index,:));
end

max_iter = 1500;
options = struct('err',1,'update',1,'test',1,'class',0,'L2',0);
tol = 10^(-12);
[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options,test_X,test_y);
fit = trop_nvar_polyval(X,num_coeffs,d) - trop_nvar_polyval(X,den_coeffs,d);


disp('L2 error')
disp(norm(train_data - fit,2))
disp('L infinity error')
disp(norm(train_data - fit,Inf))


fig15 = figure(15);
semilogy(1:max_iter,prob_out.err,'k')
hold on 
semilogy(1:max_iter,prob_out.update,'r')
semilogy(1:max_iter,prob_out.test,'b')
legend('Training error','Update norm','Testing error')
title('Error in degree 3 approximation of 6 variable function, N = 10,000')

set(fig15,'Units','Inches');
pos = get(fig15,'Position');
set(fig15,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig15,'6var_alternate_convergence_with_test','-dpdf','-r0')
hold off
%matlab2tikz('tikz/6var_alternate_convergence_with_test.tex','showInfo',false)


eval_on_test = trop_nvar_polyval(test_X,num_coeffs,d) - trop_nvar_polyval(test_X,den_coeffs,d);
%disp('L2 error')
%disp(norm(test_y-eval_on_test,2))
disp('Test L infinity error')
disp(norm(test_y - eval_on_test,Inf))

%% 10-variable problem

nvar = 10;
dmax = 1;

rng(34586) %for reproducibility 10-var

f = @(x) x(1)*x(2)*x(3) + 2*x(4)*x(5)^2*sin(x(6)) - exp(x(7)*x(8)*x(9)*x(10));
d = dmax*ones(nvar,1);

N = 10000; X = rand(N,nvar); 

train_data = zeros(N,1);
for index = 1:N
    train_data(index) = f(X(index,:));
end

test_X = rand(N,nvar); 

test_y = zeros(N,1);
for index = 1:N
   test_y(index) = f(test_X(index,:));
end

max_iter = 1500;
options = struct('err',1,'update',1,'test',1,'class',0,'L2',0);
tol = 10^(-12);

[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options,test_X,test_y);
fit = trop_nvar_polyval(X,num_coeffs,d) - trop_nvar_polyval(X,den_coeffs,d);

disp('L2 error')
disp(norm(train_data - fit,2))
disp('L infinity error')
disp(norm(train_data - fit,Inf))

fig16 = figure(16);
semilogy(1:prob_out.iterations,prob_out.err(1:prob_out.iterations),'k')
hold on 
semilogy(1:prob_out.iterations,prob_out.update(1:prob_out.iterations),'r')
semilogy(1:prob_out.iterations,prob_out.test(1:prob_out.iterations),'b')
legend('Training error','Update norm')%,'Testing error')
title('Error in degree 1 approximation of 10 variable function, N = 10,000')

%set(fig16,'Units','Inches');
%pos = get(fig16,'Position');
%set(fig16,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(fig16,'10var_alternate_convergence','-dpdf','-r0')
hold off
matlab2tikz('tikz/10var_alternate_convergence_with_test.tex','showInfo',false)

eval_on_test = trop_nvar_polyval(test_X,num_coeffs,d) - trop_nvar_polyval(test_X,den_coeffs,d);
%disp('L2 error')
%disp(norm(test_y-eval_on_test,2))
disp('Test L infinity error')
disp(norm(test_y - eval_on_test,Inf))

%% Section 9: Recovery of Univariate and Bivariate Tropical Rational Functions

Train_err = zeros(15,2); Test_err = zeros(15,2);
ntrials = 10;

for nvar = 1:2
    for dmax = 1:15
        for trial = 1:ntrials
            d = dmax*ones(nvar,1);
            
            %random tropical rational function
            true_num_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
            true_den_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
            true_num_coeffs(1) = 1;
            
            %Training Data
            N = 10000; X = 10*rand(N,nvar) - 5*ones(N,nvar); 
            train_data = trop_nvar_polyval(X,true_num_coeffs,d) - trop_nvar_polyval(X,true_den_coeffs,d);
            
            %uncomment to add noise
            %train_data = train_data + 0.05*randn(size(train_data));
            
            %Testing Data
            test_X = 10*rand(N,nvar) - 5*ones(N,nvar); 
            test_y = trop_nvar_polyval(test_X,true_num_coeffs,d) - trop_nvar_polyval(test_X,true_den_coeffs,d);
                        
            max_iter = 1000;
            options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
            tol = 10^(-12);
            
            [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options);
            fit = trop_nvar_polyval(X,num_coeffs,d) - trop_nvar_polyval(X,den_coeffs,d);                        
            eval_on_test = trop_nvar_polyval(test_X,num_coeffs,d) - trop_nvar_polyval(test_X,den_coeffs,d);
            
            Train_err(dmax,nvar) = Train_err(dmax,nvar)+(1/ntrials)*norm(train_data - fit,Inf);
            Test_err(dmax,nvar) = Test_err(dmax,nvar)+(1/ntrials)*norm(test_y - eval_on_test,Inf);
        end
    end
end

Train_err
Test_err

writematrix(Train_err,'Recovery_Degrees_Train.csv')
writematrix(Test_err,'Recovery_Degrees_Test.csv')

%% Plots for Section 9

Train_plot = readmatrix("Recovery_Degrees_Train.csv");
Test_plot = readmatrix('Recovery_Degrees_Test.csv');

fig20 = figure(20);
plot(Train_plot(:,1),'r')
hold on 
plot(Train_plot(:,2),'k')
legend('Univariate','Bivariate')
title('Training Error')
ylabel('Error')
xlabel('Degree')

set(fig20,'Units','Inches');
pos = get(fig20,'Position');
set(fig20,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig20,'low_var_train_recover','-dpdf','-r0')
hold off

fig21 = figure(21);
plot(Test_plot(:,1),'r')
hold on
plot(Test_plot(:,2),'k')
legend('Univariate','Bivariate')
title('Testing Error')
ylabel('Error')
xlabel('Degree')

set(fig21,'Units','Inches');
pos = get(fig21,'Position');
set(fig21,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig21,'low_var_test_recover','-dpdf','-r0')
hold off

%% Section 10: Recovery of 6 variable Tropical Rational Functions

Train_err6 = zeros(5,1); Test_err6 = zeros(5,1);
nvar = 6;ntrials = 5;


for dmax = 1:5
    for trial = 1:ntrials
        disp(dmax)
        disp(trial)
        d = dmax*ones(nvar,1);
                
        %random tropical rational function
        true_num_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
        true_den_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
        true_num_coeffs(1) = 1;

        
        %Training Data
        N = 10000; X = 10*rand(N,nvar) - 5*ones(N,nvar); 
        train_data = trop_nvar_polyval(X,true_num_coeffs,d) - trop_nvar_polyval(X,true_den_coeffs,d);
        
        %uncomment to add noise
        %train_data = train_data + 0.05*randn(size(train_data));
        
        %Testing Data
        test_X = 10*rand(N,nvar) - 5*ones(N,nvar); 
        test_y = trop_nvar_polyval(test_X,true_num_coeffs,d) - trop_nvar_polyval(test_X,true_den_coeffs,d);
                
        max_iter = 1000;
        options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
        tol = 10^(-8);
            
        [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options);
        fit = trop_nvar_polyval(X,num_coeffs,d) - trop_nvar_polyval(X,den_coeffs,d);
        disp('iterations:')
        disp(prob_out.iterations)
        
        %disp('L infinity error')
        %disp(norm(train_data - fit,Inf))
                
        eval_on_test = trop_nvar_polyval(test_X,num_coeffs,d) - trop_nvar_polyval(test_X,den_coeffs,d);
        %disp('Test L infinity error')
        %disp(norm(test_y - eval_on_test,Inf))
        
        Train_err6(dmax) = Train_err6(dmax)+(1/ntrials)*norm(train_data - fit,Inf)/norm(train_data,Inf);
        Test_err6(dmax) = Test_err6(dmax)+(1/ntrials)*norm(test_y - eval_on_test,Inf)/norm(test_y,Inf);
    end
end

Train_err6
Test_err6

%writematrix(Train_err6,'Recovery_Degrees_Train6.csv')
%writematrix(Test_err6,'Recovery_Degrees_Test6.csv')

%% Section 11: Recovery of 6 variable Tropical Rational Functions (with noise)

Train_err6 = zeros(5,1); Test_err6 = zeros(5,1);
nvar = 6;ntrials = 5;


for dmax = 1:2
    for trial = 1:ntrials
        disp(dmax)
        disp(trial)
        d = dmax*ones(nvar,1);
                
        %random tropical rational function
        true_num_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
        true_den_coeffs = 10*rand((dmax + 1)^nvar,1) - 5;
        true_num_coeffs(1) = 1;

        
        %Training Data
        N = 10000; X = 10*rand(N,nvar) - 5*ones(N,nvar); 
        train_data = trop_nvar_polyval(X,true_num_coeffs,d) - trop_nvar_polyval(X,true_den_coeffs,d);
        
        %uncomment to add noise
        train_data = train_data + 0.05*randn(size(train_data));
        
        %Testing Data
        test_X = 10*rand(N,nvar) - 5*ones(N,nvar); 
        test_y = trop_nvar_polyval(test_X,true_num_coeffs,d) - trop_nvar_polyval(test_X,true_den_coeffs,d);
                
        max_iter = 1000;
        options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
        tol = 10^(-8);
            
        [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options);
        fit = trop_nvar_polyval(X,num_coeffs,d) - trop_nvar_polyval(X,den_coeffs,d);
        disp('iterations:')
        disp(prob_out.iterations)

        disp('recovery difference numerator:')
        disp(num_coeffs - true_num_coeffs)
        %display('true num coeffs:')
        %display(true_num_coeffs)

        disp('recovery difference denominator:')
        disp(den_coeffs-true_den_coeffs)
        %display('true den coeffs:')
        %display(true_den_coeffs)
        
        %disp('L infinity error')
        %disp(norm(train_data - fit,Inf))
                
        eval_on_test = trop_nvar_polyval(test_X,num_coeffs,d) - trop_nvar_polyval(test_X,den_coeffs,d);
        %disp('Test L infinity error')
        %disp(norm(test_y - eval_on_test,Inf))
        
        Train_err6(dmax) = Train_err6(dmax)+(1/ntrials)*norm(train_data - fit,Inf)/norm(train_data,Inf);
        Test_err6(dmax) = Test_err6(dmax)+(1/ntrials)*norm(test_y - eval_on_test,Inf)/norm(test_y,Inf);
    end
end

Train_err6
Test_err6

%writematrix(Train_err6,'Recovery_Degrees_Train6.csv')
%writematrix(Test_err6,'Recovery_Degrees_Test6.csv')

%% Section 12: More Bivariate Examples

N = 500;

parameters_rat = zeros(6,1);
errors_rat_avg = zeros(6,1);
errors_rat_worst = zeros(6,1);
ntrials = 25;

for dmax = 1:6
for trial = 1:ntrials
X = 2*rand(N,2) -ones(N,2);
train_data = X(:,1).^2 + X(:,2).^2 + 0.25*randn(N,1);

[XX,YY] = meshgrid(-1:0.01:1,-1:0.01:1);
XY = [reshape(XX,[numel(XX),1]) reshape(YY,[numel(YY),1])];
%ZZ = XX.^2 + YY.^2; 

%K = 100;
%dmax = floor(sqrt(K/2)) - 1;


parameters_rat(dmax) = 2*(dmax+1)^2;

d = [dmax,dmax];
max_iter = 10000;
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);

%Find and Evaluate Tropical Rational Approximation
[~, ~,prob_out] = trop_nvar_rat_fit(X,train_data,max_iter,d,tol,options);
errors_rat_avg(dmax) = errors_rat_avg(dmax) + (1/ntrials)*prob_out.err(prob_out.iterations);
errors_rat_worst(dmax) = max(errors_rat_worst(dmax), prob_out.err(prob_out.iterations));
end
end


%Plot comparison with error from Maragos & Theodosis 2020
parameters_poly = [10; 25; 50; 100];
errors_poly = [0.8011; 0.6369; 0.6534; 0.6414];

figure(22)
semilogx(parameters_poly,errors_poly,'ko-')
hold on
semilogx(parameters_rat,errors_rat_avg,'bx-')
semilogx(parameters_rat,errors_rat_worst,'rs-')
legend('MT20','Rational Approximation Avg', 'Rational Approximation Worst')
xlabel('Number of Parameters')
ylabel('\|error\|_\infty')
hold off

%write plot data to csv
PlotDataPoly = table(parameters_poly,errors_poly);
PlotDataRat = table(parameters_rat,errors_rat_avg,errors_rat_worst);
writetable(PlotDataPoly,'BivariateConvexPoly.csv')
writetable(PlotDataRat,'BivariateConvexRat.csv')


%uncomment to plot the data and approximation
% figure(22)
% plot3(X(:,1),X(:,2),train_data,'k.');
% hold on
% fit = trop_nvar_polyval(XY,num_coeffs,d) - trop_nvar_polyval(XY,den_coeffs,d);
% surf(XX,YY,reshape(fit,[size(XX,1), size(YY,2)]));
% hold off
% %legend('Approximation','data')
% %title('Degree (5,5) Approximation of max(5x+ y+1,2x + 3y+3,2), noise')
% %matlab2tikz('RationalFunctionExample.tex','showInfo',false)

%% Examples from Kazda & Li

%% x^2 - y^2 example
x_lb = 0.5; x_ub = 7.5;
y_lb = 0.5; y_ub = 3.5;

N_vec = [9^2, 19^2, 49^2, 99^2];

dmax = 25;
error_by_degree = zeros(dmax,1);
time_by_degree = zeros(dmax,1);

[XX_test,YY_test] = meshgrid(x_lb:0.01:x_ub,y_lb:0.01:y_ub);
ZZ_test = XX_test.^2 - YY_test.^2;

test_x_vec = reshape(XX_test,[numel(XX_test), 1]);
test_y_vec = reshape(YY_test,[numel(YY_test), 1]);


for i = 1:4
    N = N_vec(i);
    
%     x_points = x_lb + (x_ub - x_lb)*rand(N,1);
%     y_points = y_lb + (y_ub - y_lb)*rand(N,1);
%     
%     z = exp(-10*(x_points.^2 - y_points.^2).^2);
%     data = [x_points y_points];
    
    [XX,YY] = meshgrid(x_lb:((x_ub - x_lb)/sqrt(N)):x_ub,y_lb:((y_ub - y_lb)/sqrt(N)):y_ub);
    ZZ = XX.^2 - YY.^2;
    
    x_vec = reshape(XX,[numel(XX), 1]);
    y_vec = reshape(YY,[numel(YY), 1]);
    

    for dtrial = 1:dmax
        dtrial
        d = [dtrial,dtrial];
        max_iter = 10000;
        options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
        tol = 10^(-12);
        
        tic
        [num_coeffs,den_coeffs,~] = trop_nvar_rat_fit([x_vec y_vec],reshape(ZZ,[numel(ZZ),1]),max_iter,d,tol,options);
        time_by_degree(dtrial,i) = toc;
        
        eval_on_test = trop_nvar_polyval([test_x_vec test_y_vec],num_coeffs,d) - trop_nvar_polyval([test_x_vec test_y_vec],den_coeffs,d);
        error_by_degree(dtrial,i) = norm(reshape(ZZ_test,[numel(ZZ_test),1]) - eval_on_test,'inf');
        
    end
end

figure(25)
semilogy(1:dmax, error_by_degree)
hold on
ylabel('Error (log scale)')
xlabel('Degree')
title('Error and Degree')
legend('N=100','N = 400', 'N=2500', 'N = 10000','Location','NorthEast')
hold off

figure(26)
semilogy(1:dmax,time_by_degree)
hold on 
ylabel('Time (s)')
xlabel('Degree')
title('Time and Degree')
legend('N=100','N = 400', 'N=2500', 'N = 10000','Location','NorthWest')
hold off

%Plotting Data
degree = 1:dmax;
degree = degree';
error100 = error_by_degree(:,1);
error400 = error_by_degree(:,2);
error2500 = error_by_degree(:,3);
error10000 = error_by_degree(:,4);

time100 = time_by_degree(:,1);
time400 = time_by_degree(:,2);
time2500 = time_by_degree(:,3);
time10000 = time_by_degree(:,4);

PlotData = table(degree,error100,error400,error2500,error10000,time100,time400,time2500,time10000);
writetable(PlotData,'BivariateDC.csv')

%% y^2 sin(x)/x
x_lb = 1; x_ub = 3;
y_lb = 1; y_ub = 2;

N_vec = [9^2, 19^2, 49^2, 99^2];

dmax = 25;
error_by_degree = zeros(dmax,1);
time_by_degree = zeros(dmax,1);

[XX_test,YY_test] = meshgrid(x_lb:0.01:x_ub,y_lb:0.01:y_ub);
ZZ_test = (YY_test.^2).*sin(XX_test)./XX_test;

test_x_vec = reshape(XX_test,[numel(XX_test), 1]);
test_y_vec = reshape(YY_test,[numel(YY_test), 1]);


for i = 1:4
    N = N_vec(i);
    
%     x_points = x_lb + (x_ub - x_lb)*rand(N,1);
%     y_points = y_lb + (y_ub - y_lb)*rand(N,1);
%     
%     z = exp(-10*(x_points.^2 - y_points.^2).^2);
%     data = [x_points y_points];
    
    [XX,YY] = meshgrid(x_lb:((x_ub - x_lb)/sqrt(N)):x_ub,y_lb:((y_ub - y_lb)/sqrt(N)):y_ub);
    ZZ = (YY.^2).*sin(XX)./XX;
    
    x_vec = reshape(XX,[numel(XX), 1]);
    y_vec = reshape(YY,[numel(YY), 1]);
    

    for dtrial = 1:dmax
        dtrial
        d = [dtrial,dtrial];
        max_iter = 10000;
        options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
        tol = 10^(-12);
        
        tic
        [num_coeffs,den_coeffs,prob_out] = trop_nvar_rat_fit([x_vec y_vec],reshape(ZZ,[numel(ZZ),1]),max_iter,d,tol,options);
        time_by_degree(dtrial,i) = toc;
        
        eval_on_test = trop_nvar_polyval([test_x_vec test_y_vec],num_coeffs,d) - trop_nvar_polyval([test_x_vec test_y_vec],den_coeffs,d);
        error_by_degree(dtrial,i) = norm(reshape(ZZ_test,[numel(ZZ_test),1]) - eval_on_test,'inf');
        
    end
end

figure(25)
semilogy(1:dmax, error_by_degree)
hold on
ylabel('Error (log scale)')
xlabel('Degree')
title('Error and Degree')
legend('N = 100','N = 400', 'N = 2500', 'N = 10000','Location','NorthEast')
hold off

figure(26)
semilogy(1:dmax,time_by_degree)
hold on 
ylabel('Time (s)')
xlabel('Degree')
title('Time and Degree')
legend('N = 100','N = 400', 'N = 2500', 'N = 10000','Location','NorthWest')
hold off

%Plotting Data
degree = 1:dmax;
degree = degree';
error100 = error_by_degree(:,1);
error400 = error_by_degree(:,2);
error2500 = error_by_degree(:,3);
error10000 = error_by_degree(:,4);

time100 = time_by_degree(:,1);
time400 = time_by_degree(:,2);
time2500 = time_by_degree(:,3);
time10000 = time_by_degree(:,4);

PlotData = table(degree,error100,error400,error2500,error10000,time100,time400,time2500,time10000);
writetable(PlotData,'BivariateSin.csv')

%% Example with large ridge

x_lb = 1; x_ub = 2;
y_lb = 1; y_ub = 2;

N_vec = [9^2, 19^2, 49^2, 99^2];

dmax = 25;
error_by_degree = zeros(dmax,1);
time_by_degree = zeros(dmax,1);

for i = 1:4
    N = N_vec(i);
    
    x_points = x_lb + (x_ub - x_lb)*rand(N,1);
    y_points = y_lb + (y_ub - y_lb)*rand(N,1);
    
    z = exp(-10*(x_points.^2 - y_points.^2).^2);
    data = [x_points y_points];
    
    [XX,YY] = meshgrid(x_lb:((x_ub - x_lb)/sqrt(N)):x_ub,y_lb:((y_ub - y_lb)/sqrt(N)):y_ub);
    ZZ = exp(-10*(XX.^2 - YY.^2).^2);
    
    x_vec = reshape(XX,[numel(XX), 1]);
    y_vec = reshape(YY,[numel(YY), 1]);
    
    [XX_test,YY_test] = meshgrid(x_lb:0.01:x_ub,y_lb:0.01:y_ub);
    ZZ_test = exp(-10*(XX_test.^2 - YY_test.^2).^2);
    
    test_x_vec = reshape(XX_test,[numel(XX_test), 1]);
    test_y_vec = reshape(YY_test,[numel(YY_test), 1]);
    
    for dtrial = 1:dmax
        dtrial
        d = [dtrial,dtrial];
        max_iter = 10000;
        options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
        tol = 10^(-12);
        
        tic
        [num_coeffs,den_coeffs,prob_out] = trop_nvar_rat_fit([x_vec y_vec],reshape(ZZ,[numel(ZZ),1]),max_iter,d,tol,options);
        time_by_degree(dtrial,i) = toc;
        
        eval_on_test = trop_nvar_polyval([test_x_vec test_y_vec],num_coeffs,d) - trop_nvar_polyval([test_x_vec test_y_vec],den_coeffs,d);
        error_by_degree(dtrial,i) = norm(reshape(ZZ_test,[numel(ZZ_test),1]) - eval_on_test,'inf');
        
    end
end

figure(25)
semilogy(1:dmax, error_by_degree)
hold on
ylabel('Error (log scale)')
xlabel('Degree')
title('Error and Degree')
legend('N=100','N = 400', 'N=2500', 'N = 10000','Location','NorthEast')
hold off

figure(26)
semilogy(1:dmax,time_by_degree)
hold on 
ylabel('Time (s)')
xlabel('Degree')
title('Time and Degree')
legend('N=100','N = 400', 'N=2500', 'N = 10000','Location','NorthWest')
hold off

%Plotting Data
degree = 1:dmax;
degree = degree';
error100 = error_by_degree(:,1);
error400 = error_by_degree(:,2);
error2500 = error_by_degree(:,3);
error10000 = error_by_degree(:,4);

time100 = time_by_degree(:,1);
time400 = time_by_degree(:,2);
time2500 = time_by_degree(:,3);
time10000 = time_by_degree(:,4);

PlotData = table(degree,error100,error400,error2500,error10000,time100,time400,time2500,time10000);
writetable(PlotData,'BivariateRidge.csv')
