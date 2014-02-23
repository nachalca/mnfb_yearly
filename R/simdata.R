
# Program to create the simulated data set

library(plyr)
library(mnormt)
set.seed(1234)
# Simple scenario: Z = (Z1,Z2) ~ N(0, S) 
# mu1=mu2=0, sig1=sig2=sigma and cor(Z1,Z2)=rho. 
simdat <- function(ns,r,s=1) {
  # r: is the rho23 value
  # s: is the value of the variance
  sigma <- diag(c(s,s)); sigma[2] <- r*s ; sigma[3] <- r*s
  mu <- rep(0,2)
  data.frame(rmnorm(ns, mean=mu, varcov=sigma))
}  

prm <- expand.grid(r= c(0,.25,.5,.75,.99), ns=c(10,50,250))

#prm <- expand.grid(r= c(.1,.5), s=1, n=c(50))
simdata <-  rdply(5,mdply(prm, simdat))
colnames(simdata)[1] <- 'sim'

s=list(.1,1,10,100) 

aux <- ldply(s, function(sx) {
  data.frame(simdata[,1:3],s=sx,X1=simdata$X1*sx,X2=simdata$X2*sx)
      }
)
simdata <- aux
# number of data set: 
xx <- alply(simdata[,1:4], .margins=2, unique)
prod(laply(xx, dim)[,1])

# Run from R folder 
save(simdata, file='data/simdata.Rdata')
