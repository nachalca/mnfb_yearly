# program to save a reduced piece of output to copy from home folder
library(plyr)
library(reshape2)
library(rstan)
set_cppo(mode = "fast")

d <- list.files('../data', 'sims_')
for ( i in 1:length(d) ) load(paste('../data/', d[i],sep=''))

x <- rbind(res_size10d2[[1]],res_size50d2[[1]],res_size250d2[[1]], res_size10d10[[1]],res_size50d10[[1]])
n1 <- dim(res_size10d2[[1]])[1]
n2 <- dim(res_size10d10[[1]])[1]
reduced.res  <- data.frame(dim=c(rep(2,n1*3),rep(10,n2*2)),x)
write.table(reduced.res, file='../data/reduced_res.csv', row.names=FALSE)

# checking times ..
ms=c('iw', 'siw', 'ss', 'ht')
time <- rbind(res_size10d2[[2]],res_size50d2[[2]],res_size250d2[[2]],res_size10d10[[2]],res_size50d10[[2]])
colnames(time) <- ms
df <- data.frame(n=c(10,50,250,10,50), dim=c(2,2,2,10,10), time)
write.csv(df, file='../data/timetable.csv', row.names=FALSE)

# compute the percentile that true value represent in posterior distribution
# and also the probabilty arround the true value
moreres <- function(x,nms) {
  true.r <- nms[i,'r']
  i <<- i+1
  rs <- extract(x, pars='rho', permuted = TRUE, inc_warmup=FALSE)$rho
  data.frame(per = sum(rs <= true.r)/length(rs), 
             prob= sum( (rs<=true.r+.1) & (rs >= true.r-.1) )/length(rs), 
             mae = mean( abs(rs-true.r)) )
}

pr <- c('iw', 'siw', 'ss', 'ht')
out <- NULL
for (j in 3:6) {
   nms1 <- attributes(res_size10d2[[j]])$split_labels
   nms2 <- attributes(res_size50d2[[j]])$split_labels
   nms3 <- attributes(res_size250d2[[j]])$split_labels
   nms4 <- attributes(res_size10d10[[j]])$split_labels
   nms5 <- attributes(res_size50d10[[j]])$split_labels
    i <- 1;   df1 <- data.frame(prior=pr[j-2],dim=2, ldply(res_size10d2[[j]], moreres, nms=nms1) )
   i <- 1;   df2 <- data.frame(prior=pr[j-2],dim=2 ,ldply(res_size50d2[[j]], moreres, nms=nms2) )
   i <- 1;   df3 <- data.frame(prior=pr[j-2],dim=2 ,ldply(res_size250d2[[j]], moreres, nms=nms3) )
   i <- 1;   df4 <- data.frame(prior=pr[j-2],dim=10, ldply(res_size10d10[[j]], moreres, nms=nms4) )
   i <- 1;   df5 <- data.frame(prior=pr[j-2],dim=10, ldply(res_size50d10[[j]], moreres, nms=nms5) )
   df <- rbind(df1,df2,df3,df4,df5)
   out <- rbind(out,df)
}
write.csv(out, file='../data/moreres.csv', row.names=FALSE)

#Run the IW model with scaled data on worst scenario
load('data/models_cpp.Rdata')
load('data/simdata.Rdata')
prms <- c('s1', 's2', 'rho')

runstan.sim <- function(d, it = 1200, ch = 3, w=200, prm=NULL) {
  K <- ncol(d[,-c(1:4)])
  dat = list(y = d[,-c(1:4)] , N = nrow(d), R = diag(K), k=K)
  sampling(object=m_iw, data = dat,pars=prm, iter = it, chains = ch, warmup=w)
}
printresult <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}

simd2 <- subset(simdata.2, s==0.1)
simd10 <- subset(simdata.10, s==0.1)

 # compute mean and sd to save it
momd2  <- ddply(simdata.2, .(s,sim,r,ns), function(xx) c(dim(xx)[1],apply(xx[,-c(1:4)],2,mean),apply(xx[,-c(1:4)],2,sd) ) ) 
colnames(momd2)[-c(1:4)] <- c('n',paste(rep(c('mean','sd'),each=2), rep(1:2,2), sep='.'))
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







