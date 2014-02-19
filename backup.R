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

# Get simulated data and expanded to include prior type
# load('data\\simdata.Rdata')
load('../data/simdata.Rdata')
ms=c('iw', 'siw', 'ss', 'ht')
simdata <- data.frame( ms=rep(ms, each=nrow(simdata)),rbind(simdata,simdata,simdata,simdata) )

# Compile the model objects 
source('codemodels.R')
m_iw  <- stan_model(model_code=sim.iw)
m_siw <- stan_model(model_code=sim.siw)
m_ss  <- stan_model(model_code=sim.ss)
m_ht  <- stan_model(model_code=sim.ht)


# Run simulations 
# functions to run stan model
runstan.sim <- function(d, it = 1200, ch = 3, w=200) {
  if (d$ms[1]=='iw')  mod<- m_iw
  if (d$ms[1]=='siw') mod<- m_siw
  if (d$ms[1]=='ss')  mod<- m_ss
  if (d$ms[1]=='ht')  mod<- m_ht
  dat = list(y = d[,c('X1','X2')] , N = nrow(d), R = diag( ncol(d[,c('X1','X2')])) )
  sampling(object=mod, data = dat, iter = it, chains = ch, warmup=w)
}

printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}
# testing the simulations and fitting
#zsim <- with(simdata, simdata[r==-.85 & s==1 & n==50, ])
#qplot(data=zsim, x=X1,y=X2, size=I(.5)) + geom_density2d()
# run a stan model with very few iter for compiling c code 
#dsim <- list(y = zsim[,c('X1','X2')] , N = nrow(zsim), R = diag( ncol(zsim[,c('X1','X2')])) )
#res  <- sampling(object = m_iw, data  = dsim, iter = 10, chains = 1, warmup=0)
# grid values for simulating parameters
#N <- 5
#prm <- expand.grid(r=c(-.8, -.3, 0, .3, .8), s=c(.1, 1, 50))
#prm <- expand.grid(r=c(.8), s=c(.1, 1))

simdata <- subset(simdata, sim==1 & r==.99 & s==.1 & ns==10 )                       

ptm <- proc.time()
mod_iw <-  dlply(simdata[simdata$ms=='iw', ], .(sim,r,s,ns),runstan.sim)                        
time.iw <- proc.time() - ptm
save(mod_iw,time.iw, file='../data/simula_iw.Rdata')

ptm <- proc.time()
mod_siw <-  dlply(simdata[simdata$ms=='siw', ], .(sim,r,s,ns),runstan.sim)                        
time.siw <- proc.time() - ptm
save(mod_siw,time.siw, file='../data/simula_siw.Rdata')

ptm <- proc.time()
mod_ss <-  dlply(simdata[simdata$ms=='ss', ], .(sim,r,s,ns),runstan.sim)                        
time.ss <- proc.time() - ptm
save(mod_ss,time.ss, file='../data/simula_ss.Rdata')

ptm <- proc.time()
mod_ht <-  dlply(simdata[simdata$ms=='ht', ], .(sim,r,s,ns),runstan.sim)                        
time.ht <- proc.time() - ptm
save(mod_ht,time.ht, file='../data/simula_ht.Rdata')

#res_sim <- ldply(mod_sim, printresult)
# save the results
#save(mod_sim,res_sim,time, m_iw,m_siw,m_ss,m_ht, file='../data/simulations.Rdata')
