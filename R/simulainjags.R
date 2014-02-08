# some simulations in jags
# y~ N(xb, sigma.e), b~ N(mu, Sigma), sigma.e~ IG(alpha/2, alpha*lambda/2)
library(plyr)
library(reshape2)
library(ggplot2)
library(rjags)
library(gridExtra)
library(mnormt)

# IW model
sim.jg.iw = "
model {
  for (i in 1:N) { y[i,1:2] ~ dmnorm(mu , Tau) }
  
  # Priors.  
  Tau   ~ dwish(R, df)
  Sigma  <- inverse(Tau)
  for (i in 1:2) { mu[i] ~ dnorm(0,0.001) }
  df     <- 3  
}
"
# function to run the jags models with simulated data
runjags.sim <- function(d, mod) {
  dat = list(y = d[,c('X1','X2')] , N = nrow(d), R = diag( ncol(d[,c('X1','X2')])) )
  m = jags.model(textConnection(mod), dat, n.chains=3, n.adapt=500)
  update(m, 1000)
  coda.samples(m, 2000)
}

# get the simulated data
load('data\\simdata.Rdata')

simres <- function(n,r,s, ms) {  
  sdat <- simdata[r==r & s==s & n==n, ]
  llply(ms, function(x) runjags.sim(sdat, mod=x))
}  

testing <- simres(r=-.85 , s=1 , n=50, ms=list(iw=sim.jg.iw) )


ptm <- proc.time()
res <-  mlply(prm, simres, n=200,ms=list('iw','siw')) 
proc.time() - ptm








