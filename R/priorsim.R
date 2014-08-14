# Prior Simulation from all 4 models
library(rstan)
library(plyr)
library(reshape2)
library(ggplot2)
library(devtools)


priorcode = "
  data {
  int <lower=0> k;
  int <lower=0> l;
  matrix[k,k] R;
}
parameters {
  cov_matrix[k] Q1;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta2;
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
  real<lower=0> s3_iw;
  real rho_iw;
  real rho23_iw;
  real<lower=0> s1_siw;
  real<lower=0> s2_siw;
  real<lower=0> s3_siw;
  real rho_siw;
  real rho23_siw;
  real<lower=0> s1_ss;
  real<lower=0> s2_ss;
  real<lower=0> s3_ss;
  real rho_ss;
  real rho23_ss;
  real<lower=0> s1_ht;
  real<lower=0> s2_ht;
  real<lower=0> s3_ht;
  real rho_ht;
  real rho23_ht;

// ################  IW strategy ####################
  s1_iw  <- sqrt(Sig_iw[1,1]);
  s2_iw  <- sqrt(Sig_iw[2,2]);
  rho_iw <- Sig_iw[1,2]/(s1_iw*s2_iw);
  s3_iw  <- sqrt(Sig_iw[l,l]);
  rho23_iw <- Sig_iw[2,l]/(s2_iw*s3_iw);
//====================================

// ################  sIW strategy ####################
  s1_siw  <- sqrt( sQ[1,1]*delta[1]*delta[1]);
  s2_siw  <- sqrt(sQ[2,2]*delta[2]*delta[2]);
  rho_siw <-  (sQ[1,2]*delta[2]*delta[1]) /(s1_siw*s2_siw);
  s3_siw  <- sqrt( sQ[l,l]*delta[l]*delta[l]);
  rho23_siw <- (sQ[2,l]*delta[2]*delta[l]) /(s2_siw*s3_siw);
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
        Sig_ss[m,n] <- delta2[m] * delta2[n] * Q[m,n]; 
      }
    }
    for (n in 1:k) 
      for (m in (n+1):k) 
        Sig_ss[m,n] <- Sig_ss[n,m];
  s1_ss  <- sqrt(Sig_ss[1,1]);
  s2_ss  <- sqrt(Sig_ss[2,2]);
  s3_ss  <- sqrt(Sig_ss[l,l]);
  rho_ss <- Sig_ss[1,2]/(s1_ss*s2_ss);
  rho23_ss <- Sig_ss[2,l]/(s2_ss*s3_ss);
// =================================================
// ############## HT ########################
  for (i in 1:k)   xi[i] <- 1/xi1[i];
  D <- 4*diag_matrix(xi);
  s1_ht <- sqrt(Sig_ht[1,1]);
  s2_ht<- sqrt(Sig_ht[2,2]);
  s3_ht<- sqrt(Sig_ht[l,l]);
  rho_ht <- Sig_ht[1,2]/(s1_ht*s2_ht);
  rho23_ht <- Sig_ht[2,l]/(s2_ht*s3_ht);
// ============================================
}
model {
  real mss;
  real A2;
  A2 <- 1/(1.04*1.04);
  mss <- log(.72)/2;               // SS parameter to match median of sigmas
  Sig_iw ~ inv_wishart(k+1, R);   // IW prior 
  sQ ~ inv_wishart(k+1, 0.8*R); //   sIW prior

  Q1 ~ inv_wishart( k+1, R); // correlation for SS prior
  for ( i in 1:k) delta[i] ~ lognormal(0, 1); // prior for sd on SIW prior
  for ( i in 1:k) delta2[i] ~ lognormal(mss, 1); // prior for sd on SS prior
  for (i in 1:k)  xi1[i] ~ inv_gamma(0.5, A2) ; // HT sd priors
  Sig_ht ~ inv_wishart(k+1, D); // HT prior
}
"

m_prior <- stan_model(model_code=priorcode)

simula.prs <- function(k, obj, it) {
  l <- ifelse( k > 2, 3, 2 )
  dat = list( R = diag(k), k=k, l=l)
  pr <- c(paste(rep(c('s1', 's2','s3','rho','rho23'),4), rep(c('iw','siw','ss','ht'),each=5),sep='_' ))
  res.prior <- sampling(object=obj, data = dat, pars=pr, iter=it, warmup=40)
  prior.sim <- extract(res.prior, permute=FALSE)
  d <- prior.sim
  dim(d) <- c( prod(dim(prior.sim)[1:2])  ,  dim(prior.sim)[3])  
  colnames(d) <- attributes(prior.sim)$dimnames$parameters
  data.frame(d)
}
x <- simula.prs(k=3, obj=m_prior, it=500 ) 

res <- mdply( data.frame(k=c(2,10)), simula.prs, obj=m_prior, it=2540) 
write.csv(res, file='data/priorsim.Rdata',row.names=FALSE)
res <- read.csv('data/priorsim.Rdata',header=T)
#qplot(data=d, x=log(s1_ht), y=rho_ht)
#qplot(data=d, rho_iw, geom='density')


dm <- cbind(melt(res[,c(1,grep('s1', colnames(res))) ], id.vars='k'),
            melt(res[,c(1,grep('s2', colnames(res))) ], id.vars='k'),
            melt(res[,c(1,grep('rho_', colnames(res))) ], id.vars='k'),
            melt(res[,c(1,grep('s3', colnames(res))) ], id.vars='k'),
            melt(res[,c(1,grep('rho23', colnames(res))) ], id.vars='k') 
            )

dm$p<- factor(rep(c('IW','SIW', 'BMMmu', 'HIWht'), each= dim(dm)[1]/4))
dm$prior <- factor( dm$p , levels=levels(dm$p)[c(3,4,2,1)] )
colnames(dm) <- c('dim','x1','s1','dim2','x2','s2','dim3','x3','rho','dim4','x4','s3','dim5','x5','rho23','p','prior')

# Compare quantiles....
ddply(dm, .(dim, prior), function(x) quantile(x$s2^2, c(.025, .5, .975)))
library(stats)
library(MCMCpack)

# match moment to get parameter simulations... 
# variance distributions: 
# IW(v,l) : s1^2 ~ invXi(v-d+1, l/(v-d+1)), v=d+1, l=1
# SIW(v,l,dlt) : s1^2 = xi^2*phi, 
#                xi ~ logNorm(0, 2*dlt), phi~invXi(v-d+1, l/(v-d+1))
# HT(v,A's):  s1 ~ ht(v, A), v=df, A = scale parameter
# SS(v,l,d0, d1): si^2 ~ logNorm(2*d0, 2*d1)
log(1/qgamma( c(.975,.5,.025) , shape=1, rate=1/2),10 )
log( exp(qnorm(c(.025,.5,.975), mean=log(.72)/2, sd=1))^2 ,10 )
log(  qt(c(.5125, 0.75, 0.9875),df=2, ncp=(.72/2)^2 )^2 ,10 )

htmatch <- function(A,d) {
  (1/qgamma(c(.5) , shape=1, rate=1/2) - (A*qt(c(0.75),df=2))^2 )^2 
}
# with d=2, A=1.040212
optimize(htmatch, interval=c(0, 5), d=2)$minimum
(1.04*qt(c(0.75) ,df=2))^2

# SIW is the hardest one ... 
# dirty optimization in order to match median
siwmatch <- function(pr,d=2,n=100000) {
  l <- pr[1] ; dlt <- pr[2]
  #q <- 1/qgamma( .5 , shape=1, rate=1/2)
  xi <- exp(rnorm(n, mean=0, sd=dlt))^2
  Q  <- 1/rgamma(n, shape=1, rate=l/2)
  qt <- quantile(xi*Q, probs=.5)
  round((qt- .72)^2,2) 
  #qplot(x,f1,geom='line') + geom_line(aes(x,f2), color=I('red'))
}
siwmatch(c(1,1) )
optim(par = c(1,2),fn=siwmatch, d=2)

# match the numerator in the SIw mean (which actually dont exists for v=d+1) with the 
# IW median value 0.72 :
# exp(4b1)*l = 0.72 => 4b1 + log(l) = log(0.72)
# so if b1=1 
b1 <- 1; exp(log(0.72)-4*b1)


# some plots ... 
qplot(data=res, x=rho_iw, y=rho23_iw, color=log(s2_iw,10))+scale_color_gradient2()

ddply(dm, .(dim, prior), summarise, m1 = median(s1), m2 = median(s2), v1=var(s1), v2=var(s2)) 

qplot(data=dm, x=log(s1), geom='density',color=prior, facets=~dim)
qplot(data=dm, x=log(s1), y=log(s2), size=I(.5), color=abs(rho)) + facet_grid( facets=prior~dim, scale='free')
qplot(data=dm, x=log(s1,10), geom='histogram') + facet_grid(facets=dim~prior) 
qplot(data=dm, x=rho, geom='density') + facet_grid(facets=dim~prior) 


# relation among s1 and s2 and the effect on rho
dm$dim2 <- factor(dm$dim, levels=c(2,10), labels=paste(c(2,10),'d',sep='-'))

pdf('report/figs/prior_sis2.pdf', height=5)
p <- qplot(data=dm, x=s2, y=s1,size=I(.8),color=abs(rho)) + facet_grid(facets=dim2~prior,labeller=label_parsed) 
p + xlab( expression( sigma[1] )  ) + ylab(  expression( sigma[2] )) + scale_x_log10() + scale_y_log10(limits=c(0.01,100)) +
scale_colour_gradient(name=expression('|'~ rho ~ '|'), low='white',high='red', limits=c(0,1)) +
theme(legend.position='bottom')
dev.off()

# Prior plot: countour?? how many lines ?? smooth is better ?? 
pdf('report/figs/priorsim2d.pdf', height=5)
p <- qplot(data=dm, x=s1, y=rho,size=I(.6)) + facet_grid(facets=dim2~prior, scales='free', labeller=label_parsed) 
p + xlab( expression( sigma[1] )  ) + ylab(  expression( rho )) + scale_x_log10()
dev.off()        


# gg <- ggplot_build(p)
# gg1 <- gg
# gr <- paste( rep(1:8,each=4), c('003','002', '001'), sep='-')
# gg$data[[2]]<-subset(gg1$data[[2]],group %in% gr)
# p1<-ggplot_gtable(gg)
# plot(p1)

# correlation scatter plot
qplot(data=subset(res,prior='iw' &  )
dm <- ddply(dm, .(dim,prior), mutate, ss = cut(log(s2,10), quantile(log(s2,10),probs=c(0,1/3,2/3,1)) ,include.lowest=TRUE, labels=1:3))
dm$ss <- as.numeric(as.character(dm$ss))
pdf('report/figs/priorsroro.pdf', height=4)
qplot(data=subset(dm, dim==10), x=rho, y=rho23,size=I(.8),color=ss) + facet_grid(facets=.~prior)+ 
  scale_color_gradient2(name=expression(sigma[2]),low='black',high='red', limits=c(1,3),midpoint=2, breaks=c(1,2,3),labels=c('low','med','high')) +
  theme(legend.position='bottom') + xlab( expression( rho[12] ) ) + ylab(  expression( rho[23])  ) + scale_x_continuous(breaks=c(0,1))
dev.off()        


# for discussion, inverse gamma plot
library(MCMCpack)
x <- seq(0.0001, 1/2,,1001)
ig <- dinvgamma(x, 1,1/2)
pdf(file='report/figs/ig.pdf', height=3.5)
qplot(x,ig, geom='line',ylab='density',xlab='') + geom_vline(xintercept=0.01, color=I('red')) + theme_bw()
dev.off()

n <- 10000
sims <- NULL
for (i in 1:n) sims <- rbind(sims, data.frame( t(as.numeric(rwish(3,diag(2))))))

hist(sims$X2, prob=T, 1000)

qplot(x=sims$X2, geom='histogram', binwidth=.1)
qplot(x=sims$X4, geom='histogram', binwidth=.1)













# ==========================
slopes <- ddply(dm, .(dim,prior), function(x) {
  xx <- subset(x, rho >= 0 )
  xx$ls <- log(xx$s1)
  m0 <- lm(log(rho) ~ log(ls), data=xx)
  data.frame(slope=coef(m0)[2], pval=summary(m0)$coefficients[2,4])
  }
  )
xt <- xtable(slopes, cpation='BoxCox transform and slopes')

print(xt, file='report/boxcox.tex')
qplot(data=xx, x=log(s1,10), y=rho )
p + facet_grid(facets=dim~prior) + geom_smooth(size=(2)) 
qnt <- ddply(dm, .(prior), summarise, q90 = quantile(s1, 0.9), q10=quantile(s1, 0.1))
dm2 <- merge(dm, qnt)
dm2$sdlev <- (dm2$s1 > dm2$q10) + (dm2$s1 > dm2$q90)
dm2$sdlev <- as.factor(dm2$sdlev) 
levels(dm2$sdlev) =c('low', 'med', 'high')
qplot(data=subset(dm2, sdlev!='med'),x=rho,geom='histogram') +facet_wrap(facets = prior~sdlev, ncol=2) 

            