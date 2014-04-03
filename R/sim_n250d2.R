# just run here n=250 dim=2 to get it faster


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
load('../data/models_cpp.Rdata')

# set up the parallel for plyr functions
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(6)
}

# functions to run stan model
runstan.sim <- function(d, it = 1500,ch=3, w=500, prm=NULL) {
  if (d$ms[1]=='iw')  mod<- m_iw
  if (d$ms[1]=='siw') mod<- m_siw
  if (d$ms[1]=='ss')  mod<- m_ss
  if (d$ms[1]=='ht')  mod<- m_ht
  K <- ncol(d[,-c(1:5)])
  dat = list(y = d[,-c(1:5)] , N = nrow(d), R = diag(K), k=K, mu0 = rep(0,K))
  out <- sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
  x  <- printresult(out)
  gd <- max(x$Rhat, na.rm=T); 
  #neff <- min(x$n_eff)
  j <- 1;
  it2 <- it + 1500
  w2 <- w + 500
  neff <- min(x$n_eff[x$param!='lp__'])
  
  while ( (gd > 1.1 | neff < 500) & j < 6 ) {
    out <- sampling(object=mod, data = dat,pars=prm, iter = it2, chains = ch, warmup=w2)
    x  <- printresult(out)
    gd <- max(x$Rhat, na.rm=T)
    neff <- min(x$n_eff[x$param!='lp__'])
    print(gd)
    #neff <- min(x$n_eff)
    j <- j + 1
    it2 <- it2 + 1500
    w2 <- 1000    
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
  mod_iw <-  dlply(simdata[simdata$ms =='iw', ], .(sim,r,s,ns),runstan.sim, prm=prms, .parallel=parallel)                        
  time.iw <- proc.time() - ptm
  #save(mod_iw,time.iw, file='../data/simula_iw.Rdata')
  ptm <- proc.time()
  mod_siw <-  dlply(simdata[simdata$ms=='siw', ], .(sim,r,s,ns),runstan.sim, prm=prms, .parallel=parallel)                        
  time.siw <- proc.time() - ptm
  #save(mod_siw,time.siw, file='../data/simula_siw.Rdata')
  ptm <- proc.time()
  mod_ss <-  dlply(simdata[simdata$ms=='ss', ], .(sim,r,s,ns),runstan.sim, prm=prms, .parallel=parallel)                        
  time.ss <- proc.time() - ptm
  #save(mod_ss,time.ss, file='../data/simula_ss.Rdata')
  ptm <- proc.time()
  mod_ht <-  dlply(simdata[simdata$ms=='ht', ], .(sim,r,s,ns),runstan.sim, prm=prms, .parallel=parallel)                        
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

#Run simulations for Bivariate case
data2 <- data.frame( ms=rep(ms, each=nrow(simdata.2)),rbind(simdata.2,simdata.2,simdata.2,simdata.2) )

#res_size10d2 <- simula(size=10, data=data2)
#save(res_size10d2, file='../data/sims_n10_d2.Rdata')
#res_size50d2 <- simula(size=50, data=data2)
#save(res_size50d2, file='../data/sims_n50_d2.Rdata')

res_size250d2_backup <- simula(size=250, data=data2)
save(res_size250d2_backup, file='../data/sims_n250_d2_backup.Rdata')
remove(data2)
