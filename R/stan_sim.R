# Simulation analysis

library(plyr)
library(reshape2)
library(ggplot2)
#library(rjags)
library(mnormt)

# install rstan
# options(repos = c(getOption("repos"), rstan = "http://wiki.rstan-repo.googlecode.com/git/"))
# install.packages('rstan', type = 'source')
library(rstan)
set_cppo(mode = "fast")

# Get simulated data and expanded to include prior type
# load('data/simdata.Rdata')
load('../data/simdata.Rdata')
ms=c('iw', 'siw', 'ss', 'ht')
simdata.all <- data.frame( ms=rep(ms, each=nrow(simdata)),rbind(simdata,simdata,simdata,simdata) )

# Compile the model objects 
source('codemodels.R')
m_iw  <- stan_model(model_code=sim.iw)
m_siw <- stan_model(model_code=sim.siw)
m_ss  <- stan_model(model_code=sim.ss)
m_ht  <- stan_model(model_code=sim.ht)
save(m_iw, m_siw, m_ss, m_ht, file='../data/models_cpp.Rdata')

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
simula <- function(size) {
  simdata <- subset(simdata.all, r==.75 & sim==1 & s==1 & ns==size)                       
  ptm <- proc.time()
  mod_iw <-  dlply(simdata[simdata$ms=='iw', ], .(sim,r,s,ns),runstan.sim)                        
  time.iw <- proc.time() - ptm
  #save(mod_iw,time.iw, file='../data/simula_iw.Rdata')
  ptm <- proc.time()
  mod_siw <-  dlply(simdata[simdata$ms=='siw', ], .(sim,r,s,ns),runstan.sim)                        
  time.siw <- proc.time() - ptm
  #save(mod_siw,time.siw, file='../data/simula_siw.Rdata')
  ptm <- proc.time()
  mod_ss <-  dlply(simdata[simdata$ms=='ss', ], .(sim,r,s,ns),runstan.sim,it=2000)                        
  time.ss <- proc.time() - ptm
  #save(mod_ss,time.ss, file='../data/simula_ss.Rdata')
  ptm <- proc.time()
  mod_ht <-  dlply(simdata[simdata$ms=='ht', ], .(sim,r,s,ns),runstan.sim,it=2000)                        
  time.ht <- proc.time() - ptm
  #save(mod_ht,time.ht, file='../data/simula_ht.Rdata')

  time <- c(time.iw[3],time.siw[3],time.ss[3],time.ht[3])

  res.df <- rbind( data.frame(prior='iw',ldply(mod_iw, printresult)),
                     data.frame(prior='siw',ldply(mod_siw, printresult)),
                     data.frame(prior='ss',ldply(mod_ss, printresult)),
                     data.frame(prior='ht',ldply(mod_ht, printresult)) )
list(res.df,time,mod_iw,mod_siw,mod_ss,mod_ht)
}

# Run simulations 
res_size10 <- simula(size=10)
save(res_size10, file='data/simulations10.Rdata')

res_size50 <- simula(size=50)
save(res_size50, file='data/simulations50.Rdata')

res_size250 <- simula(size=250)
save(res_size250, file='data/simulations250.Rdata')


