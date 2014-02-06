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

# Simple scenario: Z = (Z1,Z2) ~ N(0, S) 
# mu1=mu2=0, sig1=sig2=sigma and cor(Z1,Z2)=rho. 
simdat <- function(n,r,s) {
  # r: is the rho23 value
  # s: is the value of the variance
  sigma <- diag(c(s,s)); sigma[2] <- r*s ; sigma[3] <- r*s
  mu <- rep(0,2)
  data.frame(rmnorm(n, mean=mu, varcov=sigma))
}  

runstan.sim <- function(d, mod) {
  # d: are the simulated data set
  # mod: stan object with an empty comiled model
  if (mod=='siw') model = res.siw 
  if (mod=='iw')  model  = res.iw 
  dat = list(y = d , 
             N = nrow(d), 
             ns = ncol(d),
             R = diag( ncol(d)) )
  stan(fit=model, data = dat, iter = 1200, chains = 3, warmup=200)
  #stan(fit=model, data = dat,pars=pr,iter = 100, chains = 1, warmup=0)
}
# this function simulate a set of data and fit all models in ms list to it
simres <- function(n,r,s, ms) {  
  sdat <- simdat(n,r,s)
  llply(ms, function(x) runstan.sim(sdat, x))
}  

# testing the simulations and fitting
zsim <- simdat(100, .8, 50)
qplot(data=zsim, x=X1,y=X2, size=I(.5)) + geom_density2d()

# run a stan model with very few iter for compiling c code 
dsim = list(y = zsim, N = nrow(zsim), ns = ncol(zsim), R = diag( ncol(zsim)) )
res.iw  <- stan(model_code=simple.iw, data  = dsim, iter = 10, chains = 1, warmup=0)
res.siw <- stan(model_code=simple.siw, data = dsim, iter = 10, chains = 1, warmup=0)

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













#===============================================================
#===============================================================
# y~ N(xb, sigma.e), b~ N(mu, Sigma), sigma.e~ IG(alpha/2, alpha*lambda/2)
library(plyr)
library(reshape2)
library(ggplot2)
library(rjags)
library(gridExtra)
library(mnormt)

# function to run the jags models with simulated data
runjags.sim <- function(d, mod) {
  # d: are the simulated data set
  # mod: is the 'alias' for the model character object
  
  if (mod=='siw') model = mtext.hh.siw 
  if (mod=='iw') model = simple.mod 
  dat = list(y = d$y , 
             abbrev = as.numeric(d$group) ,
             year= d$x, 
             n = nrow(d), 
             ns = nlevels(factor(d$group)),
             R = diag(3))
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=100)
  update(m, 500)
  coda.samples(m, c('alpha','lambda','rho23','rho12','rho13','mu',"sigma.be",'beta','sigma.e'), 2000)
  #coda.samples(m, 'rho23', 3000)
}

# function to simulated a data set with a similar structure of the 
# bird data. The only 'free' parameter is the correlation between 
# slope and quad term. 
simdat <- function(r) {
  # r: is the rho23 value
  
  # set parameter values: mostly based on the estimations from data
  rho23 <- r
  s <- diag(c(2, .5, .5)); s[6] <- rho23*.5^2 ; s[8] <- rho23*.5^2
  mu <- rep(0,3)
  alpha <- 2 
  lambda <- 0.1
  G <- 70
  
  # get simulated data
  beta <- rmnorm(G, varcov=s)
  sigma.e <- 1/sqrt(rgamma(G,shape=alpha/2, rate=alpha*lambda/2))
  simpar <- data.frame(group=1:G,beta, sigma.e)
  colnames(simpar) <- c('group','b1', 'b2', 'b3', 'sigma.e')
  simy_f <- function(d) {
    x <- -10:10 
    eps <- rnorm(length(x), mean=0, sd=sqrt(d$sigma.e))
    y   <- with(d, b1+b2*x+b3*x^2+eps)
    data.frame(y,x)
  }
  ddply(simpar, .(group), simy_f)
}  
  
# this function simulate a set of data and fit all models in ms list to it
simres <- function(r, ms) {  
  sdat <- simdat(r)
  fitmod <- function(mod) {
    res <- runjags.sim(sdat, mod)  
    prnam <- attributes(res[[1]])$dimnames[[2]]
    aux1 <- summary(res[, c('alpha','lambda',prnam[grep('mu', prnam)], prnam[grep('sigma.be', prnam)], 'rho23', prnam[grep('beta', prnam)])])$quantile
    out <- as.data.frame(aux1)
    colnames(out) <- paste('Q', c(2.5,25,50,75,97.5),sep='')
    s <- grep('sigma.be', prnam)[c(4,6,7)]
    aux2 <- gelman.diag(res[, prnam[-s] ])
    out$gd=aux2$mpsrf    
    out$param <- rownames(out)
    out$model <- mod
    return(out)
  }
  ldply(ms, fitmod)
}  
  

# 5 values of rho, 2 types of sigma prior, N reps of each
setwd('~\\GitHub\\mnfb_yearly\\R')
source('codemodels.R')

N <- 5
rhos <- data.frame(r=c(-.8, -.3, 0, .3, .8))
ptm <- proc.time()
sim.res <- rdply(.n=N, mdply(rhos, simres, ms=list('iw','siw')) )
proc.time() - ptm

setwd('~\\GitHub\\mnfb_yearly')
save(sim.res, file='simulations.Rdata')

#--------------------------------
# to work with one example .... 
sdat <- simdat(-.8)
res <- runjags.sim(sdat, 'iw')
prnam <- attributes(res[[1]])$dimnames[[2]]

s <- grep('sigma.be', prnam)[-c(1,5,9)]
aux2 <- gelman.diag(res[, prnam[-s] ])
plot(res[, 'sigma.be[3,3]']
  
     
     
     





