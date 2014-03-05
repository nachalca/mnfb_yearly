# program to save a reduced piece of output to copy from home folder
library(plyr)
library(reshape2)
library(rstan)
set_cppo(mode = "fast")

d <- list.files('../data', 'sims_')
for ( i in 1:length(d) ) load(d[i])

x <- rbind(res_size10d2[[1]],res_size50d2[[1]],res_size250d2[[1]],res_size10d10[[1]],res_size50d10[[1]],res_size250d10[[1]])
reduced.res  <- data.frame(dim=rep(c(2,10),each=60),subset(x, param %in% c('mu[1]', 'mu[2]', 's1', 's2', 'rho') ) )

write.table(reduced.res, file='../data/reduced_res.csv', row.names=FALSE)


# checking times ..
ms=c('iw', 'siw', 'ss', 'ht')
time <- rbind(res_size10d2[[2]],res_size50d2[[2]],res_size250d2[[2]],res_size10d10[[2]],res_size50d10[[2]],res_size250d10[[2]])
colnames(time) <- ms
df <- data.frame(n=rep(c(10,50,250),2), dim=rep(c(2,10),each=3), time)
write.csv(df, file='../data/timetable.csv', row.names=FALSE)

# Run the IW model with scaled data on worst scenario

load('../data/models_cpp.Rdata')
load('../data/simdata.Rdata')
prms <- c('mu', 's1', 's2', 'rho')

runstan.sim <- function(d, it = 1200, ch = 3, w=200, prm=NULL) {
  K <- ncol(d[,-c(1:5)])
  dat = list(y = d[,-c(1:5)] , N = nrow(d), R = diag(K), k=K)
  sampling(object=m_iw, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
}
printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}

simd2 <- subset(simdata.2, s==0.1)
simd10 <- subset(simdata.10, s==0.1)

# compute mean and sd to save it
momd2  <- ddply(simd2, .(s,sim,r,ns), function(xx) c(apply(xx[,-c(1:4)],2,sd),apply(xx[,-c(1:4)],2,mean)) ) 
colnames(momd2)[-c(1:4)] <- paste(rep(c('mean','sd'),each=2), rep(1:2,2), sep='.')
momd10  <- ddply(simd10, .(s,sim,r,ns), function(xx) c(apply(xx[,-c(1:4)],2,sd),apply(xx[,-c(1:4)],2,mean)) ) 
colnames(momd10)[-c(1:4)] <- paste(rep(c('mean','sd'),each=10), rep(1:10,2), sep='.')

# scale each data set
simd2  <- ddply(simd2, .(s,sim,r,ns) , function(xx) scale(xx[,-c(1:4)]))
simd10 <- ddply(simd10, .(s,sim,r,ns) , function(xx) scale(xx[,-c(1:4)]))

mod_sciwd2  <-  dlply(simd2, .(sim,r,ns), runstan.sim, prm=prms)
mod_sciwd10 <-  dlply(simd10, .(sim,r,ns), runstan.sim, prm=prms)

res.scIW <- rbind( data.frame(dim=2,ldply(mod_sciwd2, printresult)),
                   data.frame(dim=10,ldply(mod_sciwd10, printresult)) )
                                               
save(momd2, momd10, file='../data/momentIW.Rdata')
save(mod_sciwd10, mod_sciwd2, file='../data/mod_scIW.Rdata')
write.csv(res.scIW, file='../data/res_scIW.csv', row.names=FALSE)







