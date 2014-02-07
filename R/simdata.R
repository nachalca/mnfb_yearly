
# Program to create the simulated data set

library(plyr)
library(gridExtra)
library(mnormt)

# Simple scenario: Z = (Z1,Z2) ~ N(0, S) 
# mu1=mu2=0, sig1=sig2=sigma and cor(Z1,Z2)=rho. 
simdat <- function(n,r,s) {
  # r: is the rho23 value
  # s: is the value of the variance
  sigma <- diag(c(s,s)); sigma[2] <- r*s ; sigma[3] <- r*s
  mu <- rep(0,2)
  data.frame(rmnorm(n, mean=mu, varcov=sigma))
}  

#prm <- expand.grid(r= c(0,seq(-.95, .95, .1)), s=c(.1, 1,25,50,100), n=c(50,200,500))

prm <- expand.grid(r= c(.1,.5), s=1, n=c(50))
simdata <-  mdply(prm, simdat)
save(simdata, file='simdata.Rdata')