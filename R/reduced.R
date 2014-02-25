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