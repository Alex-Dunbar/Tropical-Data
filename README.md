# Tropical-Data
MATLAB Files for Tropical Algebra and Data Project. The files are divided into two folders, Functions and Experiments.

## Functions

This folder contains the MATLAB functions used to compute tropical rational approximations to data.

Contents:

bivar_trop_polyfit.m -- Fit infinity norm minimizing bivariate tropical polynomial to given data

trop_bivar_polyval.m -- Evaluate a bivariate tropical polynomial

trop_bivar_rat_fit.m -- Use an alternating minimization method to approximate given data with a bivariate tropical rational function

trop_mat_dual_mult.m -- compute min-plus matrix vector multiplication

trop_mat_mult.m      -- compute max-plus matrix vector multiplication

trop_poly_subfit.m   -- find optimal subsolution for univariate polynomial fitting

trop_polyval.m       -- Evaluate a univariate tropical polynomial

trop_rat_fit.m       -- Use an alternating minimization method to approximate given data with a univariate tropical rational function

Still needs work:

trop_rat_2_fit.m     -- Fit a tropical rational function to data using an approximate gradient descent on the squared l2 loss

trop_rat_2_loss_grad.m -- Approximate the gradient of the squared l2 loss of a tropical rational approximation to data.


## Experiments

This folder contains the MATLAB functions used in Experiments for fitting tropical rational functions to data. 

Contents:

Bivariate_Trop_Regression_Experiments.m -- Experiments with bivariate tropical rational regression

PolyfitExperiments.m                    -- Reproduce tropical polynomial fitting examples in Maragos and Theodosis 2019 Tropical Geometry and Piecewise Linear Approximation on Weighted Lattices

RatFitExperiments.m                     -- Experiments with univariate tropical rational regression
