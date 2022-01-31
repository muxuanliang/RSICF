# R Package: RSICF
# Estimation and Inference for the Relative Contrast Function under a Monotonic Single-index Model Assumption

Version: 1.0.0

Author: Muxuan Liang <mliang@fredhutch.org>

Maintainer: Muxuan Liang <mliang@fredhutch.org>

Description: The package implements a semiparametric efficient inference of the relative contrast function under a monotonic single-index model assumption. User can choose bespoke methods to estimate outcome model and propensity. The package also provides a variaty of methods to estimate these nuisance parameters. It returns the coefficeint estimates and its sd for constructing CI of the single-index.

License: MIT

Imports: 
         glmnet,
         VariableScreening,
         mgcv,
         nloptr,
         spline2

Encoding: UTF-8

# Details
This package (function rsfit) implements the inference for coefficients of interest in a relative contrast function (see reference).

# Example

```
### Generate data

n <- 500
p <- 4
x <- matrix(runif(n*p, -0.5,0.5), c(n,p))
beta_inter <- c(1,-1,1,-1)
trt <- apply(x, 1, function(t){rbinom(1, 1, prob = exp(0.2*(t[1]^2+t[2]^2+t[1]*t[2]))/(1+exp(0.2*(t[1]^2+t[2]^2+t[1]*t[2]))))})
e <- 0.1*rnorm(n,1)
inter_effect <- (pnorm(x %*% beta_inter)-0.5)
main_effect <- sqrt(apply(x,1,function(t){sum(t^2)}))
y <- (trt-0.5) * inter_effect*main_effect + main_effect+e

### Fit our approach

fit <- rsfit(x, y, trt, splitIndex = NULL, propensityModel = 'kernel', outcomeModel = 'kernel', lossType = 'logistic', parallel = FALSE, constraint = TRUE, boundaryPoint = c(-8,8), tol = 1e-3, efficient = TRUE, local = FALSE)

### Prediction

xtest <- matrix(runif(10^5*p, -0.5,0.5), c(10^5,p)) %*% cor_matrix
    predict_trt <- as.numeric(predict.rsfit(fit$fit, newx = xtest)>0)
    
```
# References
Muxuan Liang, Menggang Yu (2022). Relative Contrast Estimation and Inference for Treatment Recommendation.

