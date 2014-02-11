
# Program to create the simulated data set

library(plyr)
library(gridExtra)
library(mnormt)
set.seed(1234)
# Simple scenario: Z = (Z1,Z2) ~ N(0, S) 
# mu1=mu2=0, sig1=sig2=sigma and cor(Z1,Z2)=rho. 
simdat <- function(ns,r,s) {
  # r: is the rho23 value
  # s: is the value of the variance
  sigma <- diag(c(s,s)); sigma[2] <- r*s ; sigma[3] <- r*s
  mu <- rep(0,2)
  data.frame(rmnorm(ns, mean=mu, varcov=sigma))
}  

prm <- expand.grid(r= round(c(-seq(.05,1,.2),0,seq(.05,1,.2)),2), s=c(.1, 1,10,50,100), ns=c(50,500))

#prm <- expand.grid(r= c(.1,.5), s=1, n=c(50))
simdata <-  rdply(5,mdply(prm, simdat))
colnames(simdata)[1] <- 'sim'
# Run from R folder 
save(simdata, file='../data/simdata.Rdata')