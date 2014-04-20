# Sampling from SIW prior, study marginal distribution to study conditional 
# distribution of rho | s1,s2

library(plyr)
library(MCMCpack)
library(ggplot2)
# set up the parallel for plyr functions
parallel <- require(doMC, quietly=TRUE)
if(parallel){
  registerDoMC(8)
}

sim.siw <- function(s1,s2,n=100) {
  d = 2
  e = 1.01
  rdply(n, {
    repeat {
      Q  <-  rwish(d+1,diag(nrow=d))
      xi <- exp(rnorm(2))
      W <- diag(xi)%*% Q %*% diag(xi)
      if (all(s1/e < diag(W)[1], diag(W)[1] < s1*e,s2/e < diag(W)[2], diag(W)[2] < s2*e)) break;
    }
    data.frame(rho=W[1,2]/sqrt(W[1,1]*W[2,2]))}
  )  
}
sig <- expand.grid(s1=c(.1,1,10),s2=c(.1,1,10))
rho.cndsiw <- mdply(sig, sim.siw, n=5000, .parallel=parallel)

save(rho.cndsiw, file = '../data/rho_cndsiw.Rdata')




