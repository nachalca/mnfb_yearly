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
  # parameters of interest 
  s1 <- sqrt(Sigma[1,1])
  s2 <- sqrt(Sigma[2,2])
  rho <- Sigma[1,2]/(s1*s2)
}
"
# functions to run the jags models with simulated data
runjags.sim <- function(d, mod, prs) {
  dat = list(y = d[,c('X1','X2')] , N = nrow(d), R = diag( ncol(d[,c('X1','X2')])) )
  m = jags.model(textConnection(mod), dat, n.chains=3, n.adapt=500)
  update(m, 1000)
  coda.samples(m, prs, 2000)
}

simres <- function(n,r,s, ms,prs) {
  # ms argument is the name of the model as character
  sdat <- simdata[r==r & s==s & n==n, ]
  if (ms=='iw') mod <- sim.jg.iw
  runjags.sim(sdat, mod, prs)
}  

# summarize models results
getresults <- function(reslist, prs) {
  #res <- unlist(reslist)
  qff <- function(x) {
    q <- summary(x[,prs])$quantile
    g <- gelman.diag(x[, prs])$psrf[,2]
    out <- data.frame(q,g)
    colnames(out) <- c(paste('Q', c(2.5,25,50,75,97.5),sep=''),'gd.upp')
    return(out)
  }
  ldply(reslist,qff )
}

# get the simulated data, and the grid for running models
load('../data/simdata.Rdata')
#load('data\\simdata.Rdata')

# testing ...
#testing  <- simres(r=-.85 , s=1 , n=50, ms='iw', prs=pars)
#prm <- expand.grid(n=50,r=c(-85, .85), s=1, ms='iw')
#testing<-  mlply(prm, simres) 
#save(testing, file='testing.Rdata')
#resutest <- getresults(testing, prs=pars)

# for real ...
prm <- expand.grid(n=unique(simdata$n),r=unique(simdata$r), s=unique(simdata$s), ms='iw')
pars <- c('mu[1]','mu[2]', 's1','s2','rho','Sigma[1,1]','Sigma[2,2]','Sigma[1,2]')

ptm <- proc.time()
models.iw <-  mlply(prm, simres) 
res.iw <- getresults(testing, prs=pars)
time.iw <- proc.time() - ptm








