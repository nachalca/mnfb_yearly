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

# Compile the model objects 
source('codemodels.R')
m_iw  <- stan_model(model_code=sim.iw)
m_siw <- stan_model(model_code=sim.siw)
m_ss  <- stan_model(model_code=sim.ss)
m_ht  <- stan_model(model_code=sim.ht)
save(m_iw, m_siw, m_ss, m_ht, file='../data/models_cpp.Rdata')

# functions to run stan model
runstan.sim <- function(d, it = 1500, ch = 3, w=500, prm=NULL) {
  if (d$ms[1]=='iw')  mod<- m_iw
  if (d$ms[1]=='siw') mod<- m_siw
  if (d$ms[1]=='ss')  mod<- m_ss
  if (d$ms[1]=='ht')  mod<- m_ht
  K <- ncol(d[,-c(1:5)])
  dat = list(y = d[,-c(1:5)] , N = nrow(d), R = diag(K), k=K, mu0 = rep(0,K))
  out <- sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
  x  <- printresult(out)
  gd <- max(x$Rhat); j <- 1;
  it <- 1500
  
  while (gd > 1.1 & j < 6) {
    it <- it + 1500
    out <- sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
    x  <- printresult(out)
    gd <- max(x$Rhat)
    j <- j + 1
  }
  return(out)
}
printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}
getiter <- function(xx) {
  attributes(xx)$stan_args[[2]]$iter
}
simula <- function(size, data) {
  prms <- c('s1', 's2', 'rho')
  simdata <- subset(data, ns==size)                       
  ptm <- proc.time()
  mod_iw <-  dlply(simdata[simdata$ms =='iw', ], .(sim,r,s,ns),runstan.sim, prm=prms)                        
  time.iw <- proc.time() - ptm
  #save(mod_iw,time.iw, file='../data/simula_iw.Rdata')
  ptm <- proc.time()
  mod_siw <-  dlply(simdata[simdata$ms=='siw', ], .(sim,r,s,ns),runstan.sim, prm=prms)                        
  time.siw <- proc.time() - ptm
  #save(mod_siw,time.siw, file='../data/simula_siw.Rdata')
  ptm <- proc.time()
  mod_ss <-  dlply(simdata[simdata$ms=='ss', ], .(sim,r,s,ns),runstan.sim, prm=prms)                        
  time.ss <- proc.time() - ptm
  #save(mod_ss,time.ss, file='../data/simula_ss.Rdata')
  ptm <- proc.time()
  mod_ht <-  dlply(simdata[simdata$ms=='ht', ], .(sim,r,s,ns),runstan.sim, prm=prms)                        
  time.ht <- proc.time() - ptm
  #save(mod_ht,time.ht, file='../data/simula_ht.Rdata')

  time <- c(time.iw[3],time.siw[3],time.ss[3],time.ht[3])

  res.df <- rbind( data.frame(prior='iw',ldply(mod_iw, printresult)),
                     data.frame(prior='siw',ldply(mod_siw, printresult)),
                     data.frame(prior='ss',ldply(mod_ss, printresult)),
                     data.frame(prior='ht',ldply(mod_ht, printresult)) )
  
  diag.df <- rbind( data.frame(prior='iw',ldply(mod_iw, getiter)),
                   data.frame(prior='siw',ldply(mod_siw, getiter)),
                   data.frame(prior='ss',ldply(mod_ss, getiter)),
                   data.frame(prior='ht',ldply(mod_ht, getiter)) )
  
  out <- list(res.df,time,mod_iw,mod_siw,mod_ht,mod_ss,diag.df)
  names(out) <- c('res', 'times', 'iw', 'siw', 'ht', 'ss','diag')
  return(out)
}

# Get simulated data and expanded to include prior type
# load('data/simdata.Rdata')
load('../data/simdata.Rdata')
ms=c('iw', 'siw', 'ht', 'ss')

# Run simulations for Bivariate case
# data2 <- data.frame( ms=rep(ms, each=nrow(simdata.2)),rbind(simdata.2,simdata.2,simdata.2,simdata.2) )
# res_size10d2 <- simula(size=10, data=data2)
# save(res_size10d2, file='../data/sims_n10_d2.Rdata')
# res_size50d2 <- simula(size=50, data=data2)
# save(res_size50d2, file='../data/sims_n50_d2.Rdata')
# res_size250d2 <- simula(size=250, data=data2)
# save(res_size250d2, file='../data/sims_n250_d2.Rdata')
# remove(data2)

# Run simulations for 10 dimension case: takes too long, so reduce the simulations
# only size: 10 and 50
# only rho: 0, .99
# only s: .1, 1, 100

d <- data.frame( ms=rep(ms, each=nrow(simdata.10)),rbind(simdata.10,simdata.10,simdata.10,simdata.10) )
data10 <- subset(d, s %in% c(.1,1,100) & r %in% c(0,.99))

res_size10d10 <- simula(size=10, data=data10)
save(res_size10d10, file='../data/sims_n10_d10.Rdata')

res_size50d10 <- simula(size=50, data=data10)
save(res_size50d10, file='../data/sims_n50_d10.Rdata')

#res_size250d10 <- simula(size=250, data=data10)
#save(res_size250d10, file='../data/sims_n250_d10.Rdata')
remove(data10)


# Run simulations for 10 dimension case
#data100 <- data.frame( ms=rep(ms, each=nrow(simdata.100)),rbind(simdata.100,simdata.100,simdata.100,simdata.100) )
#res_size10 <- simula(size=10, data=data100)
#save(res_size10, file='../data/sims_n10_d100.Rdata')
#res_size50 <- simula(size=50, data=data100)
#save(res_size50, file='../data/sims_n50_d100.Rdata')
#res_size250 <- simula(size=250, data=data100)
#save(res_size250, file='../data/sims_n250_d100.Rdata')
#remove(data100)