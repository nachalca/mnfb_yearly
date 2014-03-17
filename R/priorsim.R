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
  matrix[k,k] D1;
  matrix[k,k] Q; 
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
  for (i in 1:k) {
    delta1[i] <- 1/sqrt(Q1[i,i]);
  }
//  for (n in 1:k) {
//    for (m in 1:(n-1) ) {
//      Q[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
//    }
//  }
//  for (n in 1:k) 
//    for (m in (n+1):k) 
//      Q[m,n] <- Q[n,m];
    D1 <- diag_matrix(delta1);
    Q <- D1 * Q1 * D1;
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

simula.prs <- function(k, obj, it) {
  dat = list( R = diag(k), k=k)
  pr <- paste( rep(c('s1', 's2', 'rho'),4), rep(c('iw','siw','ss','ht'),each=3),sep='_' )
  res.prior <- sampling(object=obj, data = dat, pars=pr, iter=it, warmup=40)
  prior.sim <- extract(res.prior, permute=FALSE)
  d <- prior.sim
  dim(d) <- c( prod(dim(prior.sim)[1:2])  ,  dim(prior.sim)[3])  
  colnames(d) <- attributes(prior.sim)$dimnames$parameters
  data.frame(d)
}

x <- simula.prs(k=100, obj=m_prior, it=200 ) 

res <- mdply( data.frame(k=c(2,10)), simula.prs, obj=m_prior, it=2400 ) 

write.csv(res, file='data/priorsim.Rdata',row.names=FALSE)
res <- read.csv('data/priorsim.Rdata', header=T)

#qplot(data=d, x=log(s1_ht), y=rho_ht)
#qplot(data=d, rho_iw, geom='density')
dm <- cbind(melt(res[,c(1,grep('s1', colnames(res))) ], id.vars='k'),melt(res[,c(1,grep('rho', colnames(res))) ], id.vars='k') )
dm$p<- factor(rep(c('iw','siw', 'ss', 'ht'), each= dim(dm)[1]/4))
dm$prior <- factor( dm$p , levels=levels(dm$p)[c(2,3,1,4)] )
colnames(dm) <- c('dim','x1','s1','dim2','x2','rho','p','prior')

pdf('report/figs/priorsim2d.pdf')
qplot(data=dm, y=s1, x=rho,size=I(.5))+geom_smooth(se=FALSE, size=I(1.5), color=I('red')) + facet_wrap(facets=dim~prior, ncol=4, scale='free_y')
dev.off()        

qplot(data=dm, x=log(s1), geom='histogram')+facet_wrap(facets=~prior)


qnt <- ddply(dm, .(prior), summarise, q90 = quantile(s1, 0.9), q10=quantile(s1, 0.1))
dm2 <- merge(dm, qnt)
dm2$sdlev <- (dm2$s1 > dm2$q10) + (dm2$s1 > dm2$q90)
dm2$sdlev <- as.factor(dm2$sdlev) 
levels(dm2$sdlev) =c('low', 'med', 'high')
qplot(data=subset(dm2, sdlev!='med'),x=rho,geom='histogram') +facet_wrap(facets = prior~sdlev, ncol=2) 

            