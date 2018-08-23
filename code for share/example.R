# Example using SS for estimating a cov matrix

# 1) libraries 
library(plyr)
library(reshape2)
library(ggplot2)
library(mnormt)

# install rstan: 
# look at the instructions in webpage to install rstan
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 2)

#2) Stan code for model: y ~ N(mu, Sigma), mu ~ normal(), Sigma ~ SS[LN, IW] 

# Sigma ~ SS[LN, IW] means Sigma = D * R * D, 
# where D is diagonal with standard deviations that are log-noarmal distributed
# R is a correlatin matrix computed from an IW distribution.

# this compile and generate c++ code for the model in ss.stan
mod <- stan_model('ss.stan')


# 3) Small simulation

R <- diag(3)
R[R==0] <- .7

S <- diag(c(1,1.2,1.8)) %*% R %*% diag(c(1,1.2,1.8))
y <- rmnorm(20, mean=rep(0,3),  varcov = S)

k = ncol(y)
dat <- list(y = y, N = nrow(y), R = diag(k), k=k, mu0 = rep(0,k))

# compute the results
out <- sampling(object=mod, data = dat )

# 4) look at posterior summaries 

# delta: standard deviation
# R correlation matrix
# Sigma covariance matrix 
# the rest are auxiliary parameters
out
