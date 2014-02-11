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
  coda.samples(m, prs, 3000)
}
# summarize models results
getresults <- function(reslist, prs) {
  #res <- unlist(reslist)
  qff <- function(x) {
    q <- summary(x[,prs])$quantile
    g <- gelman.diag(x[, prs])$psrf[,2]
    out <- data.frame(prs,q,g)
    colnames(out) <- c('param',paste('Q', c(2.5,25,50,75,97.5),sep=''),'gd.upp')
    return(out)
  }
  ldply(reslist,qff )
}

# get the simulated data, and the grid for running models
load('../data/simdata.Rdata')
#load('data\\simdata.Rdata')

# testing ...
#testing  <- simres(r=-.85 , s=1 , n=50, ms='iw', prs=pars)
#prm <- expand.grid(n=50,r=r1, s=1, ms='iw')
#testing<-  mlply(prm, simres) 
#save(testing, file='testing.Rdata')
#resutest <- getresults(testing, prs=pars)
# more testing 
#simdata$r1 <- round(simdata$r, 2)
#r1 <- round(c(-seq(.05,1,.2),0,seq(.05,1,.2)),2)
#prm.t <- expand.grid(ns=unique(simdata$n),rs=r1, ss=unique(simdata$s))
#simres.test <- function(ns,rs,ss) {
# ms argument is the name of the model as character
#  sdat <- with(simdata, simdata[r1==rs & s==ss & n==ns,])
#  dim(sdat)
#}
# simres <- function(n,r,s, ms,...) {
#   # ms argument is the name of the model as character
#   sdat <- with(simdata, simdata[r1==r & s==s & n==n, ])
#   if (ms=='iw') mod <- sim.jg.iw
#   #runjags.sim(sdat, mod,...)
#   dim(sdat)
# }  
#aux <- mdply(prm.t, simres.test)
#head(aux)

# for real ...
pars <- c('mu[1]','mu[2]', 's1','s2','rho','Sigma[1,1]','Sigma[2,2]','Sigma[1,2]')
ptm <- proc.time()
models.iw <-  dlply(simdata,.(sim,r,s), runjags.sim, prs=pars, mod=sim.jg.iw) 
res.iw <- getresults(models.iw, prs=pars)
time.iw <- proc.time() - ptm
# save in data folder
save(models.iw, res.iw, time.iw, file='../data/iw.simres.Rdata')
