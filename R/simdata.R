
# Program to create the simulated data set

library(plyr)
library(mnormt)
# Simple scenario: Z = (Z1,Z2) ~ N(0, S) mu1=mu2=0, sig1=sig2=sigma and cor(Z1,Z2)=rho. 
# More dimensions: k-vector normally distributed with zero mean, equal variance and equal 
# correlation for each pair of vars.

# We first simulate with Variance = 1 and then re-scale to obtain other data sets. 
simdat <- function(ns,r,dim) {
  # r: is the common rho value for all dimensions
  sigma <- diag(dim); l <- lower.tri(sigma); u <- upper.tri(sigma)
  sigma[l] <- r; sigma[u]<- r
  mu <- rep(0,dim)
  data.frame(rmnorm(ns, mean=mu, varcov=sigma))
}  

# Simulate 5 data sets for each combination of r,size keeping separate the dimension
prm <- expand.grid(r= c(0,.25,.5,.75,.99), ns=c(10,50,250))

set.seed(1234)
s100 <-  rdply(5,mdply(prm, simdat, dim=100))
colnames(s100)[1] <- 'sim'
s2 <- s100[,1:5]
s10 <- s100[,1:13]

# test how close is pearson to the truth
corrs <- ddply(s2, .(sim,r,ns), summarise, pearson = cor(X1,X2))
qplot(data=corrs, r,pearson, color=sim, facets=~ns) + geom_abline(1)

# rescale each data set, ss represent the new standar deviation
ss <- data.frame(s=c(.01,.1,1,10,100))
rescale    <- function(sx, dx) data.frame(dx[,1:3],dx[,-c(1:3)]*sx)
simdata.2    <- mdply(ss, rescale, dx=s2)
simdata.10   <- mdply(ss, rescale, dx=s10)
simdata.100  <- mdply(ss, rescale, dx=s100)
  
# number of data set: 
#xx <- alply(simdata.2[,1:4], .margins=2, unique)
#prod(laply(xx, dim)[,1])

# Run from R folder 
save(simdata.2, simdata.10, simdata.100, file='../data/simdata.Rdata')

