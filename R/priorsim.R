# Prior Simulation from all 4 models

library(rstan)
library(plyr)
library(reshape2)
library(ggplot2)


priorcode = "
  data {
  int <lower=0> k;
  matrix[k,k] R;
}
parameters {
  cov_matrix[k] Q1;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] xi1;
  cov_matrix[k] sQ;
  cov_matrix[k] Sig_iw;
  cov_matrix[k] Sig_ht;
}
transformed parameters {
  matrix[k,k] D;
  corr_matrix[k] Q; 
  cov_matrix[k] Sig_ss;
  vector<lower=0>[k] delta1;
  vector<lower=0>[k] xi;
// I save firs rho and firs variance for each prior
  real<lower=0> s1_iw;
  real<lower=0> s2_iw;
  real rho_iw;
  real<lower=0> s1_siw;
  real<lower=0> s2_siw;
  real rho_siw;
  real<lower=0> s1_ss;
  real<lower=0> s2_ss;
  real rho_ss;
  real<lower=0> s1_ht;
  real<lower=0> s2_ht;
  real rho_ht;
// ################  IW strategy ####################
  s1_iw  <- sqrt(Sig_iw[1,1]);
  s2_iw  <- sqrt(Sig_iw[2,2]);
  rho_iw <- Sig_iw[1,2]/(s1_iw*s2_iw);
//====================================

// ################  sIW strategy ####################
  s1_siw  <- sqrt( sQ[1,1]*delta[1]*delta[1]);
  s2_siw  <- sqrt(sQ[2,2]*delta[2]*delta[2]);
  rho_siw <-  (sQ[1,2]*delta[2]*delta[1]) /(s1_siw*s2_siw);
# ===================================

// ################  SS strategy ####################
// Q is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
// a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops
  for (i in 1:k) delta1[i] <- 1/sqrt(Q1[i,i]);
  for (n in 1:k) {
    for (m in 1:n) {
      Q[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
    }
  }
  for (n in 1:k) 
    for (m in (n+1):k) 
      Q[m,n] <- Q[n,m];
// compute covariance matrix as: Sigma = D*Q*D, where D = diag(delta) 
    for (n in 1:k) {
      for (m in 1:n) {
        Sig_ss[m,n] <- delta[m] * delta[n] * Q[m,n]; 
      }
    }
    for (n in 1:k) 
      for (m in (n+1):k) 
        Sig_ss[m,n] <- Sig_ss[n,m];
  s1_ss  <- sqrt(Sig_ss[1,1]);
  s2_ss  <- sqrt(Sig_ss[2,2]);
  rho_ss <- Sig_ss[1,2]/(s1_ss*s2_ss);
// =================================================
// ############## HT ########################
  for (i in 1:k)   xi[i] <- 1/xi1[i];
  D <- 4*diag_matrix(xi);
  s1_ht <- sqrt(Sig_ht[1,1]);
  s2_ht<- sqrt(Sig_ht[2,2]);
  rho_ht <- Sig_ht[1,2]/(s1_ht*s2_ht);
// ============================================
}
model {
  Sig_iw ~ inv_wishart(k+1, R); // IW prior 

  sQ ~ inv_wishart(k+1, R); //   sIW prior 

  Q1 ~ inv_wishart( k+1, R); // correlation for SS prior
  for ( i in 1:k) delta[i] ~ lognormal(0, 1); // prior for sd on SS, SIW priors

  for (i in 1:k)  xi1[i] ~ inv_gamma(0.5, 0.1); // HT sd priors
  Sig_ht ~ inv_wishart(k+1, D); // HT prior
}
"


m_prior <- stan_model(model_code=priorcode)
dat = list( R = diag(2), k=2)
pr <- paste( rep(c('s1', 's2', 'rho'),4), rep(c('iw','siw','ss','ht'),each=3),sep='_' )
res.prior <- sampling(object=m_prior, data = dat, pars=pr)

print(res.prior)

prior.sim <- extract(res.prior, permute=FALSE)
d <- prior.sim
dim(d) <- c(4000,13)
colnames(d) <- attributes(prior.sim)$dimnames$parameters
d <- data.frame(d)
#colnames(d) <- gsub("\\.",replacement='', colnames(d))
write.csv(d, file='data/priorsim.Rdata',row.names=FALSE)
d <- read.csv('data/priorsim.Rdata', header=T)

qplot(data=d, x=log(s1_ht), y=rho_ht)
qplot(data=d, rho_iw, geom='density')

dm <- cbind(melt(d[,c(1,4,7,10)]),melt(d[,c(2,5,8,11)]),melt(d[,c(3,6,9,12)]) )
dm$prior<- rep(c('iw','siw', 'ss', 'ht'), each=4000 )
colnames(dm) <- c('x1', 's1', 'x2', 's2', 'x3', 'rho', 'prior')
pdf('report/figs/priorsim2d.pdf')
qplot(data=dm, y=log(s1), x=rho,size=I(1))+geom_smooth(se=FALSE, size=I(1.5), color=I('red')) + facet_wrap(facets=~prior, scale='free_y')
dev.off()        

qplot(data=dm, x=log(s1), geom='histogram')+facet_wrap(facets=~prior)


qnt <- ddply(dm, .(prior), summarise, q90 = quantile(s1, 0.9), q10=quantile(s1, 0.1))
dm2 <- merge(dm, qnt)
dm2$sdlev <- (dm2$s1 > dm2$q10) + (dm2$s1 > dm2$q90)
dm2$sdlev <- as.factor(dm2$sdlev) 
levels(dm2$sdlev) =c('low', 'med', 'high')
qplot(data=subset(dm2, sdlev!='med'),x=rho,geom='histogram') +facet_wrap(facets = prior~sdlev, ncol=2) 

            