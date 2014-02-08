# Simulation analysis

library(plyr)
library(reshape2)
library(ggplot2)
library(rjags)
library(gridExtra)
library(mnormt)
# install rstan
# options(repos = c(getOption("repos"), rstan = "http://wiki.rstan-repo.googlecode.com/git/"))
# install.packages('rstan', type = 'source')
library(rstan)
set_cppo(mode = "fast")

# Get simulated data
load('data\\simdata.Rdata')

# Compile the model objects 
m_iw <- stan_model(model_code=sim.iw)
m_siw <- stan_model(model_code=sim.siw)
m_ss <- stan_model(model_code=sim.ss)
m_ht <- stan_model(model_code=sim.ht)

# Run simulations 

# functions to run each simulation
runstan.sim <- function(d, mod, iter = 1200, chains = 3, warmup=200) {
  # d: are the simulated data set
  # mod: stan object with an empty comiled model
  dat = list(y = d[,c('X1','X2')] , N = nrow(d), R = diag( ncol(d[,c('X1','X2')])) )
  sampling(object=mod, data = dat, iter = iter, chains = chains, warmup=warmup)
}
# this function simulate a set of data and fit all models in ms list to it
simres <- function(n,r,s, ms) {  
  sdat <- simdat[r==r & s==s & n==n, ]
  llply(ms, function(x) runstan.sim(sdat, mod=x))
}  

# testing the simulations and fitting
zsim <- with(simdata, simdata[r==-.85 & s==1 & n==50, ])
qplot(data=zsim, x=X1,y=X2, size=I(.5)) + geom_density2d()

# run a stan model with very few iter for compiling c code 
dsim <- list(y = zsim[,c('X1','X2')] , N = nrow(zsim), R = diag( ncol(zsim[,c('X1','X2')])) )
res  <- sampling(object = m_iw, data  = dsim, iter = 10, chains = 1, warmup=0)

# grid values for simulating parameters
#N <- 5
prm <- expand.grid(r=c(-.8, -.3, 0, .3, .8), s=c(.1, 1, 50))

#prm <- expand.grid(r=c(.8), s=c(.1, 1))

ptm <- proc.time()
res <-  mlply(prm, simres, n=200,ms=list('iw','siw')) 
proc.time() - ptm

aux <- 
# save the results
res.sim <- unlist(res)
setwd('~\\GitHub\\mnfb_yearly')
save(res.sim, file='simulations.Rdata')





