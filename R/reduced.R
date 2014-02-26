# program to save a reduced piece of output to copy from home folder
load('../data/simulations10.Rdata')
load('../data/simulations50.Rdata')
load('../data/simulations250.Rdata')

reduced.res <- rbind(
  subset(res_size10[[1]], param %in% c('mu[1]', 'mu[2]', 's1','s2','rho' )),
  subset(res_size50[[1]], param %in% c('mu[1]', 'mu[2]', 's1','s2','rho' )),
  subset(res_size250[[1]], param %in% c('mu[1]', 'mu[2]', 's1','s2','rho' ))
)

write.table(reduced.res, file='../data/reduced_res.csv', row.names=FALSE)

# compute the percentile that true value represent in posterior distribution
# and also the probabilty arround the true value
moreres <- function(x) {
  true.r <- nms[i,'r']
  i <<- i+1
  rs <- extract(x, pars='rho', permuted = TRUE, inc_warmup=FALSE)$rho
  data.frame(per=sum(rs <= true.r)/length(rs), 
             prob=sum( (rs<=true.r+.1) & (rs >= true.r-.1) )/length(rs) )
}

pr <- c('iw', 'siw', 'ss', 'ht')
out <- NULL
for (j in 3:6) {
  nms1 <- attributes(res_size10[[j]])$split_labels
  nms2 <- attributes(res_size50[[j]])$split_labels
  nms3 <- attributes(res_size250[[j]])$split_labels
  i <- 1;   df1 <- data.frame(prior=pr[j-2], ldply(res_size10[[j]], moreres) )
  i <- 1;   df2 <- data.frame(prior=pr[j-2], ldply(res_size50[[j]], moreres) )
  i <- 1;   df3 <- data.frame(prior=pr[j-2], ldply(res_size250[[j]], moreres) )
  df <- rbind(df1,df2,df3)
  out <- rbind(out,df)
  }

write.table(out, file='../data/percentile.csv', row.names=FALSE)


