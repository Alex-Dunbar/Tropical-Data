%% Single Variable tropical rational function fitting tests
% Section 1: Fitting Tropical line Using alternating method
% Section 2: Fitting Tropical rational function using alternating method
% Section 3: Fitting sine curve using alternating method
% Section 4: Fitting sine curve in various degrees

%% Section 1: Piecewise linear function
x = linspace(-1,12,200); 
y = max(x-2,3);
x = x'; y = y'; 
%y1 = y + 0.25*randn(size(y)); %gaussian noise
y1 = y + rand(size(y))-0.5*ones(size(y)); %uniform noise

max_iter = 100; d = 1;
options = struct('err',0,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);
%Find and Evaluate Tropical Rational Approximation
[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(x,y1,max_iter,d,tol,options);
fit = trop_nvar_polyval(x,num_coeffs,d) - trop_nvar_polyval(x,den_coeffs,d);

%compare with polynomial regression
[~,polyfit] = trop_nvar_polyfit(x,y,d);

figure(1)
plot(x,y1,'k.')
hold on
plot(x,fit,'b');
plot(x,polyfit,'r')
title('Approximation of a Tropical Line')
legend('Data','Rational Approximation','Polynomial Approximation','Location','best')
hold off

%Write plot data to file
Plot_Data = table(x,y1,fit,polyfit);
writetable(Plot_Data,'Tropical_Line_Approx.csv')


%% Visualization of Loss Landscape for above

loss = @(p1,p2,q2) norm(y1 - max(p1,x+p2) + max(0,x+q2),'inf');
loss_mat = zeros(numel(-2:0.1:2));

i = 1; q2 = -1000000; %q2 = den_coeffs(2)-den_coeffs(1);
p1vec = -5:0.1:5;
p2vec = -5:0.1:5;
for p1 = p1vec
    j = 1;
    for p2 = p2vec
        loss_mat(j,i) = loss(p1,p2,q2);
        j = j+1;
    end
    i = i+1;
end

figure()
[X,Y] = meshgrid(p1vec,p2vec);
surf(X,Y,loss_mat)
hold on
plot3(num_coeffs(1)-den_coeffs(1),num_coeffs(2)-den_coeffs(1),loss(num_coeffs(1)-den_coeffs(1),num_coeffs(2)-den_coeffs(1),den_coeffs(2)-den_coeffs(1)),'rx')
hold off

%% Section 2: Tropical Rational function
y = max(x+2,4) - max(x-3,1); y3 = y + 0.05*randn(size(y));

max_iter = 100; d = 2;
options = struct('err',0,'update',0,'test',0,'class',0,'L2',0);
tol = 10^(-12);
%Find and Evaluate Tropical Rational Approximation
[num_coeffs, den_coeffs,~] = trop_nvar_rat_fit(x,y3,max_iter,d,tol,options);
fit = trop_nvar_polyval(x,num_coeffs,d) - trop_nvar_polyval(x,den_coeffs,d);

figure(2)
plot(x,y3,'k.')
hold on
plot(x,fit,'b')
title('Approximation of a Tropical Rational Function')
legend('data','Approximation','Location','best')
hold off


%% Section 3: Sine
x = linspace(-1,12,200); x=x';
y = sin(x); y2 = y + 0.05*randn(size(y)); 

%uncomment to save data to .csv for training a corresponding NN.
%writematrix([0;x],"sin_x.csv")
%writematrix([0;y2],"sin_y.csv")

%Many iterations/degree 15
max_iter = 10000;
d = 15;

%compute fit
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);
%Find and Evaluate Tropical Rational Approximation
[num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(x,y2,max_iter,d,tol,options);
fit = trop_nvar_polyval(x,num_coeffs,d) - trop_nvar_polyval(x,den_coeffs,d);

%uncomment to save rational function parameters for initializing NN
% writematrix([0;num_coeffs],"15sin_num.csv")
% writematrix([0;den_coeffs],"15sin_den.csv")
% degrees = 0:d;
% writematrix([0;degrees'],"15sin_exps.csv");

fig4 = figure(4);
plot(x,y2,'k.')
hold on
plot(x,fit,'b');
plot(x,sin(x),'color',[0 0.5 0])
title('Degree 15 Tropical Rational Fit')
legend('Data','Approximation','sin x','location','best')

% set(fig4,'Units','Inches');
% pos = get(fig4,'Position');
% set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig4,'sin_1515_alternate_fit','-dpdf','-r0')

matlab2tikz('tikz/sin_1515_alternate_fit.tex','showInfo',false)
hold off

fig5 = figure(5);
semilogy(1:prob_out.iterations,prob_out.err(1:prob_out.iterations),'k')
hold on
semilogy(1:prob_out.iterations,prob_out.update(1:prob_out.iterations),'r')
title('Error for Degree 15 Fit to Sin Data')
xlabel('Iteration')
ylabel('Error (log scale)')
legend('\infty-norm error','update norm','location','best')

% set(fig5,'Units','Inches');
% pos = get(fig5,'Position');
% set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig5,'sin_1515_alternate_convergence','-dpdf','-r0')

matlab2tikz('tikz/sin_1515_alternate_convergence.tex','showInfo',false)

hold off

% Uncomment to check nondifferentiability at solution
% tol = 10^(-13); %can set tolerance for difference between 
% [satisfied] = check_opt(x,y2,num_coeffs,den_coeffs,d,tol)



%% Section 4: Dependence of error on degrees

x = linspace(-1,12,200); x=x';
y = sin(x); y2 = y + 0.05*randn(size(y)); 

maxiter = 10000;
options = struct('err',1,'update',1,'test',0,'class',0,'L2',0);
tol = 10^(-12);

for d = 1:20
    %Find and Evaluate Tropical Rational Approximation
    [num_coeffs, den_coeffs,prob_out] = trop_nvar_rat_fit(x,y2,max_iter,d,tol,options);
    fit = trop_nvar_polyval(x,num_coeffs,d) - trop_nvar_polyval(x,den_coeffs,d);
    
    errors(d) = prob_out.err(prob_out.iterations);
    iters(d) = prob_out.iterations;
end

fig6 = figure(6);
plot(1:20, errors)
hold on
title('Dependence of Error on Degree')
xlabel('Degree')
ylabel('Error')

% set(fig6,'Units','Inches');
% pos = get(fig6,'Position');
% set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig6,'degree_error_sin','-dpdf','-r0')

matlab2tikz('tikz/degree_error_sin.tex','showInfo',false)
hold off

fig7 = figure(7);
plot(1:20, iters)
hold on
title('Dependence of Iterations on Degree')
xlabel('Degree')
ylabel('Number of iterations')

% set(fig7,'Units','Inches');
% pos = get(fig7,'Position');
% set(fig7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig7,'degree_iter_sin','-dpdf','-r0')


matlab2tikz('tikz/degree_iter_sin.tex','showInfo',false)
hold off

