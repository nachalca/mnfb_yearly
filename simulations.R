# Simulation analysis
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
  if (mod=='iw') model = mtext.hh.iw 
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

###########################

sdat <- simdat(-.8)
res <- runjags.sim(sdat, 'siw')
prnam <- attributes(res[[1]])$dimnames[[2]]
s <- grep('sigma.be', prnam)[-c(1,5,9)]
aux2 <- gelman.diag(res[, prnam[-s] ])
plot(res[, 'sigma.be[3,3]'])






