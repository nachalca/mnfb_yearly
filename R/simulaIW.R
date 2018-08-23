

###############################################################################
# Simulate from Inverse Wishart for different values in the lambda matrix
library(plyr)
library(ggplot2)
library(MCMCpack)
? rwish

# s1,s2 are sd of the prior matrix parameter, s12 is the prior covariance
simIW <- function(s1=1, s2=1, rho=0) {
  s12 <- rho*s1*s2
  s <- matrix(c(s1^2,s12,s12,s2^2), ncol=2)
  r <- matrix ( riwish(v=3, S=s)[c(1,2,4)] ,ncol=3) 
  data.frame( SD1=sqrt(r[,1]),SD2=sqrt(r[,3]),cor=r[,2]/sqrt(r[,1]*r[,3])  )
}

rdply(10, simIW)
sig <- expand.grid(s1=c(0.001,0.01,0.1),s2=c(0.001,0.01,0.1), rho=0)

sims <- rdply(5000, mdply(sig, simIW) )

qplot(data=sims, x=SD1, geom='density', facets=s2~s1) + scale_x_log10()
qplot(data=sims, x=cor, geom='density') + facet_grid(facets=s2~s1,scales='free')

###############################################################################

###### 22 Agosto, 2018 #### 











