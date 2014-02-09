# some simulations in jags
# y~ N(xb, sigma.e), b~ N(mu, Sigma), sigma.e~ IG(alpha/2, alpha*lambda/2)
library(plyr)
library(reshape2)
library(ggplot2)
library(rjags)
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
  m = jags.model(textConnection(mod), dat, n.chains=3, n.adapt=50)
  update(m, 100)
  coda.samples(m,c('mu', 'Sigma'), 200)
}

# get the simulated data
load('../data/simdata.Rdata')

simres <- function(n,r,s, ms) {
  # ms argument is the name of the model as character
  sdat <- simdata[r==r & s==s & n==n, ]
  runjags.sim(sdat, mod=get(ms) )
}  
# testing  <- simres(r=-.85 , s=1 , n=50, ms=sim.jg.iw )

prm <- expand.grid(r=c(-85, .85), s=1,n=50, ms='sim.jg.iw')
testing<-  mlply(prm, simres) 
save(testing, file='testing.Rdata')


#ptm <- proc.time()
#res <-  mlply(prm, simres, n=200,ms=list('iw','siw')) 
#proc.time() - ptm








