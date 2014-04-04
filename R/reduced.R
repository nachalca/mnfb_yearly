# program to save a reduced piece of output to copy from home folder
# library(plyr)
# library(reshape2)
# library(rstan)
# set_cppo(mode = "fast")
# 
# d <- list.files('../data', 'sims_')
# for ( i in 1:length(d) ) load(paste('../data/', d[i],sep=''))
# 
# #x <- rbind(res_size10d2[[1]],res_size50d2[[1]],res_size250d2[[1]], res_size10d10[[1]],res_size50d10[[1]])
# x <- rbind(res_size10d2[[1]],res_size50d2[[1]],res_size250d2_backup[[1]], res_size10d10[[1]],res_size50d10[[1]])
# 
# n1 <- dim(res_size10d2[[1]])[1]
# n2 <- dim(res_size10d10[[1]])[1]
# 
# reduced.res  <- data.frame(dim=c(rep(2,n1*3),rep(10,n2*2)),x)
# write.table(reduced.res, file='../data/reduced_res.csv', row.names=FALSE)
# 
# # get the number of iterations
# diagres <- rbind(res_size10d2$diag,res_size50d2$diag,res_size250d2_backup$diag, res_size10d10$diag,res_size50d10$diag)
# nn1 <- dim(res_size10d2$diag)[1] 
# nn2 <- dim(res_size10d10$diag)[1]
# diagres <- data.frame(dim=c(rep(2,nn1*3),rep(10,nn2*2)),diagres)
# write.table(diagres, file='../data/diag_res.csv', row.names=FALSE)
# 
# # get the times
# ms <- c('iw', 'siw', 'ss', 'ht')
# times <- rbind(res_size10d2[[2]],res_size50d2[[2]],res_size250d2_backup[[2]],res_size10d10[[2]],res_size50d10[[2]])
# colnames(times) <- ms
# df <- data.frame(n=c(10,50,250,10,50), dim=c(2,2,2,10,10), times)
# write.csv(df, file='../data/timetable.csv', row.names=FALSE)



#Run the IW model with scaled data on worst scenario
load('../data/models_cpp.Rdata')
load('../data/simdata.Rdata')
prms <- c('s1', 's2', 'rho')
# set up the parallel for plyr functions
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(4)
}

runstan.sim <- function(d, it = 1200, ch = 3, w=200, prm=NULL) {
  K <- ncol(d[,-c(1:4)])
  dat = list(y = d[,-c(1:4)] , N = nrow(d), R = diag(K), k=K,mu0 = rep(0,K))
  sampling(object=m_iw, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
}
printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}

simd2 <- subset(simdata.2, s %in% c(0.01, 0.1) )
simd10 <- subset(simdata.10, s %in% c(0.01, 0.1) )

# compute mean and sd to save it
momd2  <- ddply(simdata.2, .(s,sim,r,ns), function(xx) c(dim(xx)[1],apply(xx[,-c(1:4)],2,mean),apply(xx[,-c(1:4)],2,sd) ) ) 
colnames(momd2)[-c(1:4)] <- c('n',paste(rep(c('mean','sd'),each=2), rep(1:2,2), sep='.'))
momd10  <- ddply(simd10, .(s,sim,r,ns), function(xx) c(apply(xx[,-c(1:4)],2,sd),apply(xx[,-c(1:4)],2,mean)) ) 
colnames(momd10)[-c(1:4)] <- paste(rep(c('mean','sd'),each=10), rep(1:10,2), sep='.') 

# scale each data set
simd2  <- ddply(simd2, .(s,sim,r,ns) , function(xx) scale(xx[,-c(1:4)], center=FALSE))
simd10 <- ddply(simd10, .(s,sim,r,ns) , function(xx) scale(xx[,-c(1:4)], center=FALSE))
prms <- c('s1','s2','rho')

mod_sciwd2  <-  dlply(simd2, .(s,sim,r,ns), runstan.sim, prm=prms, parallel=parallel)
mod_sciwd10 <-  dlply(simd10, .(s,sim,r,ns), runstan.sim, prm=prms, parallel=parallel)
res.scIW <- rbind( data.frame(dim=2,ldply(mod_sciwd2, printresult)),
                   data.frame(dim=10,ldply(mod_sciwd10, printresult)) )                                               

save(momd2, momd10, file='../data/momentIW.Rdata')
save(mod_sciwd10, mod_sciwd2, file='../data/mod_scIW.Rdata')
write.csv(res.scIW, file='../data/res_scIW.csv', row.names=FALSE)


