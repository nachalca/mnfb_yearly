# RESULTS 

# PART I: Simulations results, plots and tables
res  <- read.table('data/reduced_res.csv', header=T)
rescov  <- read.csv('data/rescov.csv', header=T)
times <- read.csv('data/timetable.csv', header=T)
scIW <- read.csv('data/res_scIW.csv', header=T)
iterations <- read.table('data/diag_res.csv', header=T)

# reorder and relabel prior factor
res$prior.old <- res$prior

res$prior <- factor(res$prior, levels=c('iw','siw', 'ht', 'ss'),labels=c('IW','SIW', 'HIWht', 'BMMmu'))

rescov$prior <- factor(rescov$.id, levels=c('iw','siw', 'ht', 'ss'),labels=c('IW','SIW', 'HIWht', 'BMMmu'))

with(res, table(prior,prior.old))

# for compare against pearson instead to true value
load('data/simdata.Rdata')
corr <- ddply(simdata.2, .(s,sim,r,ns), summarize, pearson=cor(X1,X2))
res.p <- merge(res, corr)

library(ggplot2)
library(plyr)
library(reshape2)
library(rstan)
library(xtable)

# Gelman diag plot 
pdf('report/figs/gelmandiagd2.pdf')
qplot(data=subset(res,dim==2),x=param ,ymin=1,ymax=Rhat,geom='linerange',size=I(1))+facet_grid(facets=s~prior,scale='free') + geom_hline(yintercept=1.1, colour='red') 
dev.off()

ppdf('report/figs/gelmandiagd10.pdf')
qplot(data=subset(res,dim==10),x=param ,ymin=1,ymax=Rhat,geom='linerange',size=I(1))+facet_wrap(facets=s~prior,scale='free_y') + geom_hline(yintercept=1.1, colour='red') 
dev.off()

# pdf('report/figs/diagd10_it3000.pdf')
# qplot(data=subset(res,dim==10 & s==100),x=param ,ymin=1,ymax=Rhat,geom='linerange',size=I(1))+
#   facet_wrap(facets=sim~prior,scale='free_y',nrow=5) + geom_hline(yintercept=1.1, colour='red') +
#   xlab('') + ylab('Gelman-Rubin diagnostic statistic after 3000 iterations')
# dev.off()


# Plots for rhos estimated vs true
pdf('report/figs/fig_rho_d2.pdf', height=4.5)
d <- subset(res, param=='rho' & dim==2)
d$rx <- d$r + runif(nrow(d),-.1,.1)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$s <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data=d ,x=rx, y=mean,color=prior,shape=prior,xlab='True Correlation', ylab='Correlation Posterior Mean') + 
  facet_grid(facets=ns~s,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,0.5,1))
dev.off()

pdf('report/figs/fig_rho_d10.pdf', height=4.5)
d <- subset(res, param=='rho' & dim==10)
d$rx <- d$r + runif(nrow(d),-.1,.1)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$s <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data=d ,x=rx, y=X50.,color=prior,shape=prior,xlab='True Correlation', ylab='Posterior Median') + facet_grid(facets=ns~s,scales='free', labeller=label_parsed) + 
  geom_abline(1) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,0.5,1))
dev.off()

# compare results in 2 and ten dimensions
pdf('report/figs/fig_d2d10.pdf', height=4.5)
d <- subset(res, param=='rho',select=c(1:8))
d <- dcast(d, prior+sim+r+s+ns+param ~ dim)
colnames(d)[7:8] <- c('dim2','dim10')
d <- d[!is.na(d$dim10),]
d$ns <- factor(d$ns, levels=c(10,50), labels=paste('n',c(10,50),sep='=='))
d$s <- factor(d$s, levels=c(0.01,1,100), labels=paste('sigma',c(0.01,1,100),sep='=='))
qplot(data=d ,x=dim2, y=dim10,color=prior,shape=prior,size=I(3), xlab='Bivariate model', ylab='Ten-dimension model') + 
  facet_grid(facets=ns~s, labeller=label_parsed) + geom_abline(1)+ theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,1))
dev.off()

# compare against pearson instead to true value
load('data/simdata.Rdata')
corr <- ddply(simdata.2, .(s,sim,r,ns), summarize, pearson=cor(X1,X2))
res.p <- merge(subset(res,dim==2), corr)
pdf('report/figs/fig_rho_d2.pdf', height=4.5)
d <- subset(res.p, param=='rho')
d$rx <- d$r + runif(nrow(d),-.05,.05)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$s <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data=d ,x=pearson, y=mean,color=prior,shape=prior,xlab='True Correlation', ylab='Posterior Mean') + 
  facet_grid(facets=ns~s,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,0.5,1))
dev.off()


# Inferences about variances ... 
pdf('report/figs/fig_s1_d2.pdf', height=4.5)
d <- subset(res, param=='s2' & dim==2)
d$sx <- d$s*runif(nrow(d),.5, 1.5)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$r <- factor(d$r, levels=c(0.00, 0.25, 0.50, 0.75, 0.99), labels=paste('rho',c(0, 0.25, 0.5, 0.75, 0.99),sep='=='))
qplot(data=d ,x=sx, y=X50.,color=prior,shape=prior,xlab='True Standard deviation', ylab='Standard deviation posterior Mean') + 
  facet_grid(facets=ns~r,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_log10() + scale_y_log10()
dev.off()

qplot(data=d, x=r,y=mean,color=prior,geom='jitter')+ facet_grid(facets= s~ns,scale='free_y')

qplot(data=subset(d,s==100&ns=='n==50'), x=r,y=mean,color=prior)
+ scale_x_log10() + scale_y_log10()

pdf('report/figs/fig_s1_d10.pdf', height=4.5)
d <- subset(res, param=='s1' & dim==10)
d$sx <- d$s*runif(nrow(d),.5, 1.5)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$r <- factor(d$r, levels=c(0.00, 0.25, 0.50, 0.75, 0.99), labels=paste('rho',c(0, 0.25, 0.5, 0.75, 0.99),sep='=='))
qplot(data=d ,x=sx, y=mean,color=prior,shape=prior,xlab='True Standard deviation', ylab='Standard deviation posterior Mean') + 
  facet_grid(facets=ns~r,scales='free',labeller=label_parsed) + geom_abline(1) + theme(legend.position= 'bottom') + scale_x_log10() + scale_y_log10()
dev.off()

# Inference about covariance cov12
rescov$cov <- with(rescov, r*s*s)
with(rescov,  table(cov, r))

pdf('report/figs/fig_cov_d2.pdf', height=4.5)
d <- subset(rescov, dim==2)
d$covx <- d$cov*runif(nrow(d),.5, 1.5)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='='))
d$s <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('s',c(0.01,0.1,1,10,100),sep='='))
qplot(data=d ,x=cov, y=mean,color=prior,shape=prior,xlab='True Covariance', ylab='Covariance posterior Mean') + 
  facet_wrap(facets=ns~s,scales='free') + geom_abline(1) + theme(legend.position= 'bottom' ,axis.text.x=element_blank(),axis.text.y=element_blank())
dev.off()

# work with scaled data to "fix" IW inference
pdf('report/figs/scIW.pdf', height=4.5)
d <- subset(scIW, param=='rho' & dim=='2')
# d <- subset(d, Rhat < 2) # 2 models does not converge out of 300
d$rx <- d$r + runif(nrow(d),-.05,.05)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$dim2 <- factor(d$dim, levels=c(2,10), labels=paste('dim',c(2,10),sep='='))
d$ss <- factor(d$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data=d ,x=rx,color=I('red'), y=mean,xlab='True Correlation', ylab='Posterior Mean') + facet_grid(facets=ns~ss,scales='free', labeller=label_parsed) + 
  geom_abline(1) + scale_color_discrete(name=element_blank()) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,1))
dev.off()

# for this plot we need to bring the scaling factor used in each data set
pdf('report/figs/scIWs1.pdf', height=6)
d <- subset(scIW, param=='rho' & dim==2)
#d <- subset(scIW, param=='s1' & s %in% c(0.01, 0.1))
# d <- subset(d, Rhat < 2) # 2 models does not converge out of 300
d$sx <- d$s*runif(nrow(d),.5, 1.5)
d$ns <- factor(d$ns, levels=c(10,50,250), labels=paste('n',c(10,50,250),sep='=='))
d$dim2 <- factor(d$dim, levels=c(2,10), labels=paste('dim',c(2,10),sep='='))
d$r <- factor(d$r, levels=c(0.00, 0.25, 0.50, 0.75, 0.99), labels=paste('rho',c(0, 0.25, 0.5, 0.75, 0.99),sep='=='))
qplot(data=d ,x=sx,color=dim2, y=mean,xlab='True standard deviation', ylab='Posterior Mean') + facet_grid(facets=r~ns,scales='free', labeller=label_parsed) + 
  geom_abline(1) + scale_color_discrete(name=element_blank()) + theme(legend.position= 'bottom') + scale_x_continuous(breaks=c(0,1))
dev.off()


# Iterations, times, effective samples .... 
y <- melt(times, id.vars=c('n','dim'))
colnames(y)[-2] <- c('ns', 'prior', 'time')

tot <- sum(times[,3:6])
apply(times[,3:6], 2,sum)/tot


x <- ddply( subset(res,param!='lp__') , .(dim,prior,param,ns), summarise, samples = mean(n_eff))

samptab <- xtable(dcast(x, dim ~ prior), caption='Effective sample size')
print(samptab,file= 'report/samptab.tex', caption.placement='top', include.rownames=FALSE)

z <- merge(x,y)
z$ratio <- z$time/z$samples
z$ratio.rn <- round(z$time/z$samples , 3)
xt <- xtable(dcast(z, dim+ns+param~prior, value.var='ratio'),digits=2, label='time',caption='Iteration time per effective samples (seconds)')
print(xt,file= 'report/time_ratio.tex', hline.after=c(-1,9,9,15) ,caption.placement='top', include.rownames=FALSE)

pdf('report/figs/times.pdf', height=4)
#d <-subset(z,ns==10)
#d$dimension <- factor(d$dim, levels=c(2,10)) #, labels=paste('dim',c(2,10),sep='='))
z$param2 <- factor(z$param, levels=c('rho','s1','s2'), labels=c('rho','sigma[1]','sigma[2]') )
qplot(data=z,x=ns,y=ratio,color=prior) + facet_grid(facets=~dim,labeller=label_parsed)        
      xlab='Prior distribution', ylab='Time per samples (seconds)') 
+ facet_grid(facets=~param2,labeller=label_parsed) 
dev.off()


# ==========
# CI length computation, we compute it using fisher transformation for 
# correlatio

res$mae <- with(res, abs(X50. - r) )

# for rho CI
res$low <-   with(res, log( (1+X2.5.)/(1-X2.5.) )/2 )
res$upp <-   with(res, log( (1+X97.5.)/(1-X97.5.) )/2)
res$length <- with(res, upp - low)
mae.res <- ddply(subset(res,param=='rho'), .(prior,dim,s,r,ns), summarise,mae=mean(mae),length=mean(length)) 

pdf('report/figs/cilength.pdf', height=4)
mae.res$sample.size <- factor(mae.res$ns, levels=c(10,50,250) )
mae.res$s2 <- factor(mae.res$s, levels=c(0.01,0.1,1,10,100), labels=paste('sigma',c(0.01,0.1,1,10,100),sep='=='))
qplot(data = subset(mae.res,dim==2 & ns==10), y=length,x=r,color=prior,size=I(3),shape=sample.size,ylab='Credible interval length', xlab='True Correlation') + 
  geom_line() + facet_grid(facets=~s,labeller=label_parsed) + theme(legend.position='bottom') + scale_x_continuous(breaks=c(0,0.5,1)) 
dev.off()

# for sigmas CI
res$length[res$param %in% c('s1', 's2')] <- with(subset(res,param %in% c('s1', 's2')), (X97.5.-X2.5.))
mae.res2 <- ddply(subset(res,param=='s1'), .(prior,dim,s,r,ns), summarise,mae=mean(mae),length=mean(length)) 

pdf('report/figs/cilength_s1.pdf', height=4)
mae.res2$sample.size <- factor(mae.res2$ns, levels=c(10,50,250) )
mae.res2$sx <- mae.res2$s * runif(nrow(mae.res2),.5, 1.5)
mae.res2$r <- factor(mae.res2$r, levels=c(0.00, 0.25, 0.50, 0.75, 0.99), labels=paste('rho',c(0, 0.25, 0.5, 0.75, 0.99),sep='=='))
qplot(data = subset(mae.res2,dim==2 & ns==10), y=length,x=sx,color=prior,size=I(3),shape=prior,ylab='Credible interval length', xlab='True standard deviation') + 
  geom_line() + facet_grid(facets=~r,labeller=label_parsed) + theme(legend.position='bottom') + scale_x_log10() +  scale_y_log10()
dev.off()


t <- seq(0.0001,.9999,,100)
qnorm(t, sd=10)
qplot(t, qnorm(t), geom='line') + geom_line(aes(t,qnorm(t,sd=3)),color=I('red'))


pdf('report/figs/mae2d.pdf')
qplot(data=subset(mae.res, dim==2), y=mae,x=rx,color=prior,facets=ns~s,geom='smooth',method='lm',se=FALSE)
dev.off()
pdf('report/figs/mae2d.pdf')
qplot(data=subset(mae.res, dim==10), y=mae,x=rx,color=prior,facets=ns~s,geom='smooth',method='lm',se=FALSE)
dev.off()


# ================================================
# ================================================
# PART II: Real data exapmple, modelling yearly average of bird count species. 
# Data are created on mnfb repository, with the yearly_data.R code. 

bird.yeartotal <- read.csv('data/bird_yeartotal.csv')

# compute average using the count.add (this step should be on yearly_data.R)
bird.yeartotal$ave.add <- with(bird.yeartotal, count.add/samples)

# create a centered year variable
bird.yeartotal$yearc <- bird.yeartotal$year-2000
#setwd('~\\GitHub\\mnfb_yearly')

# species info
code <- read.csv('~/Documents/mnfb/data/nrri_bird_code.csv',header=T)

# libraries 
library(xtable)
library(plyr)
library(reshape2)
library(ggplot2)
library(GGally)
library(rstan)
#-------------------------------------------------
# 0) data description
tab_f <- function(x,sp) {
  x <- x[x$abbrev %in% sp,]
  data.frame(x[,c('count','abbrev')])
}

# get the 10 more abundant species
aux <- sort(with(bird.yeartotal, tapply(count,abbrev,sum)),decreasing=T )[1:10]
spn <- names(aux)
tab <- ddply(bird.yeartotal, .(year,forestN), tab_f, sp=spn)
#qplot(data=subset(tab, abbrev=='NAWA'), x=year, y=count,geom='line',color=forestN)

# pull out year 2007 to compare with the online report
aux <- reshape(subset(tab, year==2007) , direction='wide', timevar='forestN', idvar='abbrev',drop='year')

sp <- subset(code, abbrev %in% aux$abbrev, select=c('common','abbrev'))
mostab07 <- merge(sp,aux)[,c(2,1,3:5)]
colnames(mostab07) <- c('Specie','Abbrev',levels(bird.yeartotal$forestN) )
mostab07 <- mostab07[order(mostab07[,3],decreasing=T),]

tab = xtable(mostab07, caption='Total counts on 2007 for the 10 most abundant species', label='count07')
print(tab, file="report/highcounts.tex", include.rownames=FALSE, caption.placement='top')

# Overall trend within forest 
totales <- ddply(bird.yeartotal, .(year,forestN), summarise, count=sum(count) )
totales$Forest <- totales$forestN
pdf(file='report/figs/rawtrend.pdf', height=5)
qplot(data=totales,x=year,xlab='Year',ylab='Total Bird Count',y=count,shape=Forest,color=Forest, size=I(3))+
  theme(legend.position='bottom') +geom_line()
dev.off()


#===============================================================
#1) Species correlation: We want to estimate the correlation among species using: 
#    average count and total count
#    bivariate models and ten-varaite models
#    the 4 priors for the covariance matrix

# compute correlation among all species 
d <- dcast(subset(bird.yeartotal,forestN=='Superior'), year~abbrev, value.var='ave.add' )
bigcor <- cor(d[,-1])
bigcor.m <- melt(bigcor)

m <- with(bigcor.m, tapply(value, X2, mean))
mcnt <- with(subset(bird.yeartotal,forestN=='Superior'), tapply(count.add,abbrev, mean))
aux <- merge(data.frame(nm = names(m), corr=m), data.frame(nm=names(mcnt), cnt=mcnt))
qplot(data=aux, cnt, corr)
sp <- c(names(head(sort(m),5) ),names(tail(sort(m),5) ) )

colnames(bigcor.m)[-3] <- c('X1', 'X2')
bigcor.m$X2b <- reorder(bigcor.m$X2, bigcor.m$value, mean)
bigcor.m$X1b <- factor(bigcor.m$X1, levels=levels(bigcor.m$X2b))
dd <- subset(bigcor.m,X1%in% sp & X2 %in% sp)
qplot(data=dd , X1b,X2b,geom='tile', fill=value) + scale_fill_gradient2(low='black',high='red',midpoint=0, limits = c(-1, 1))


# use the 10 species most abundant species in 2007, use only superior forest data
# we work with centered data
spcorr.dtraw <- subset(bird.yeartotal, abbrev %in% mostab07$Abbrev & forestN=='Superior',
                    select= c('year','count.add',"abbrev","ave.add","yearc","samples")) 
colnames(spcorr.dtraw)[c(2,4)] <- c("count.add.or","ave.add.or" )
aux <- ddply(spcorr.dtraw, .(abbrev), summarise, 
             count.add = scale(count.add.or, scale=FALSE), 
             ave.add = scale(ave.add.or, scale=FALSE), 
             year = year,
             count.add.or = count.add.or,             
             yearc= yearc) 
spcorr.dt <- aux[order(aux$year), ]

spnm <- subset(code, abbrev %in% spcorr.dtraw$abbrev, select=c('common','abbrev'))
spcorr.dtraw<- merge(spcorr.dtraw, spnm)

tab <- ddply(spcorr.dtraw, .(common), summarise, 
             mean.count=mean(count.add.or), 
        sd.count = sd(count.add.or), 
        mean.ave =mean(ave.add.or), sd.ave = sd(ave.add.or))
tab <- tab[order(tab$mean.count,decreasing=T),]


xt <- xtable(tab, caption='Summary statistics of bird count data in Superior National Forest from
1995 to 2013 for the 10 most abundant species.', label='tab:bird', digits=c(NA,NA,0,0,2,2))
print(xt, file="report/bird_stat.tex", include.rownames=FALSE, caption.placement='top')

pdf(file='report/figs/rawtrend.pdf', height=4)
qplot(data=spcorr.dtraw,x=year, xlab='Year',ylab='Total Bird Count',y=count.add.or,color=common, size=I(3))+
  geom_line() + theme(legend.title=element_blank(), legend.position='bottom') + 
  guides(col = guide_legend(ncol = 5, byrow = TRUE))
dev.off()

pdf(file='report/figs/rawtrend_ave.pdf', height=4)
qplot(data=spcorr.dtraw,x=year, xlab='Year',ylab='Average Bird Count',y=ave.add.or,color=common, size=I(3))+
  geom_line() + theme(legend.title=element_blank(), legend.position='bottom') + 
  guides(col = guide_legend(ncol = 5, byrow = TRUE))
dev.off()


# total samples points per year
samples <- ddply(spcorr.dtraw, .(year), summarise, samples=mean(samples))
pdf(file='report/figs/samples.pdf', height=4)
qplot(data=samples,x=year, xlab='Year',ylab='Total sampling points',y=samples, size=I(3))+
  geom_line() + ylim(c(400,600))
#+ theme(legend.title=element_blank(), legend.position='bottom') + 
#  guides(col = guide_legend(ncol = 5, byrow = TRUE))
dev.off()



# We write 3 functions: f1 to create data set, f2 to fit the stan models, f3 to collect results. 
f1 <- function(dt, dim, vv) { 
  sp <- unique(dt$abbrev)
  ns <- length(sp)
  if (dim == 2) {
    # create a data set for estimating all pairwise correlation with bivariate models
    fn2d <- function(sp1, sp2, d=dt, v=vv) {
      d <- d[order(d$year), ]
      sp1 <- as.character(sp1); sp2 <- as.character(sp2)
      if (v=='ave') {
        x <- with(d, ave.add[abbrev==sp1])
        y <- with(d, ave.add[abbrev==sp2])
      }
      if (v=='count') {
        x <- with(d, count.add[abbrev==sp1])
        y <- with(d, count.add[abbrev==sp2])
      }
      df <- data.frame(sp1,sp2, year=unique(d$year), x, y)
      colnames(df) <- c( 'sp1', 'sp2','year',paste(v,1:2, sep='.'))
      return(df)
    }
    sppairs <- expand.grid(sp1=sp,sp2=sp )
    x <- matrix(1:100, ncol=10)
    sppairs <- sppairs[ x[lower.tri(x)],]
    out <- mdply(sppairs,  fn2d, d=dt)  
  }
  if (dim==10) {
    # create data set for estimating correlation on 10 data set. 
    if (vv=='ave')  out <- dcast(dt[,c('year','abbrev', 'ave.add')], year~abbrev)
    if (vv=='count') out <- dcast(dt[,c('year','abbrev', 'count.add')], year~abbrev)
    colnames(out)[-1] <- paste(vv,colnames(out)[-1], sep='.' )
    }  
  return(out)
}

f2 <- function(dt,dim,it = 1200, ch=3, w=200, prm=NULL,ms,v,scl=FALSE) {
  if (ms=='iw')  mod<- m_iw; if (ms=='siw') mod<- m_siw; if (ms=='ss')  mod<- m_ss; if (ms=='ht')  mod<- m_ht  
  if (ms == 'iw.sc')  {
    mod<- m_iw
    dt[,grep(v, colnames(dt))] <- scale(dt[,grep(v, colnames(dt))], center=FALSE)
    }
  dat = list(y = dt[,grep(v, colnames(dt))] , N = nrow(dt), R = diag(dim), k=dim, mu0=rep(0,dim))
  sampling(object=mod, data = dat,pars=prm, iter = it, chains = ch, warmup=w)    
}
f3 <- function(xx) {
  x <- data.frame(summary(xx)$summary)
  data.frame(param=rownames(x), round(x[,1:8],4),n_eff=round(x$n_eff),Rhat=x[,10])
}

# Create 4 data sets: dim=2,10 and response=average and count. 
allpairs.ave <- f1(dt=spcorr.dt, dim=2, vv='ave')
allpairs.cnt <- f1(dt=spcorr.dt, dim=2, vv='count')
all10.ave <- f1(dt=spcorr.dt, dim=10, vv='ave')
all10.cnt <- f1(dt=spcorr.dt, dim=10, vv='count')

# Run the models
load('data/models_cpp.Rdata')
mss <- c('iw', 'siw', 'ht', 'ss','iw.sc')
res.2dave  <- list(); res.2dcnt  <- list()
res.10dave <- list(); res.10dcnt <- list()
for (i in 1:5) {
  m <- mss[i]
  res.2dave[[i]] <- dlply(allpairs.ave,.(sp1,sp2),f2, dim=2,v='ave',ms=m, prm= c('s1', 's2', 'rho'))
  res.2dcnt[[i]] <- dlply(allpairs.cnt,.(sp1,sp2),f2, dim=2,v='count',ms=m, prm= c('s1', 's2', 'rho'))

  res.10dave[[i]] <- f2(dt=all10.ave, dim=10,v='ave', ms=m, prm=c('Sigma','Rho'))
  res.10dcnt[[i]] <- f2(dt=all10.cnt, dim=10,v='count', ms=m, prm=c('Sigma','Rho'))
}
names(res.10dave) <- mss ; names(res.10dcnt) <- mss
names(res.2dave) <- mss ; names(res.2dcnt) <- mss
save(res.2dave, res.2dcnt, res.10dave, res.10dcnt, file='data/res_birdcorr.Rdata')

load('data/res_birdcorr.Rdata')
# Collect results
x2dave <- ldply(res.2dave, function(x) ldply(x, f3) )
x2dcnt <- ldply(res.2dcnt, function(x) ldply(x, f3) )
x10dave <- ldply(res.10dave, f3) 
x10dcnt <- ldply(res.10dcnt, f3) 

# put names into the results from 10 dimension data
spnm <- gsub('ave.', replacement='', colnames(all10.ave)[-1])
par = outer(spnm, spnm, paste, sep='.')
dim(par) <- c(100,1)
x10dave.2 <- data.frame(pair=rep(par,5),subset(x10dave,param!='lp__' ) )
x10dcnt.2 <- data.frame(pair=rep(par,5),subset(x10dcnt,param!='lp__' ) )
  
# subset only the lower triangular of coefficients
x <- matrix(1:100, ncol=10)
pos2 <- x[lower.tri(x)]
pos1 <- diag(x)
x10dave.2 <- ddply(x10dave.2, .(.id), function(xx) xx[c(pos1,100+pos2),])
x10dcnt.2 <- ddply(x10dcnt.2, .(.id), function(xx) xx[c(pos1,100+pos2),])

# for pairwise models subset only rhos and add pair variable
x2dave$pair <- with(x2dave, paste(sp1,sp2,sep='.'))
x2dcnt$pair <- with(x2dcnt, paste(sp1,sp2,sep='.'))

# put all correlatoins toghether in one data sset
res_corr <- rbind(data.frame(type='pairwise',var='ave',x2dave[x2dave$param=='rho',c('.id', 'pair', 'mean')]),
                  data.frame(type='pairwise',var='cnt',x2dcnt[x2dcnt$param=='rho',c('.id', 'pair', 'mean')]),
                  data.frame(type='tendim',var='ave',x10dave.2[grep('Rho', x10dave.2$param),c('.id', 'pair', 'mean')]),
                  data.frame(type='tendim',var='cnt',x10dcnt.2[grep('Rho', x10dcnt.2$param),c('.id', 'pair', 'mean')])
            )
                
res_var <- rbind(
  data.frame(type='pairwise',var='ave',x2dave[x2dave$param=='s1',c('.id', 'pair', 'mean')]),
  data.frame(type='pairwise',var='cnt',x2dcnt[x2dcnt$param=='rho',c('.id', 'pair', 'mean')]),
  data.frame(type='tendim',var='ave',x10dave.2[grep('Sigma', x10dave.2$param),c('.id', 'pair', 'mean')]),
  data.frame(type='tendim',var='cnt',x10dcnt.2[grep('Sigma', x10dcnt.2$param),c('.id', 'pair', 'mean')])
)

# plot species trends and compute sample correlation and sd for each pair
sample.ave.corr <- ddply(allpairs.ave, .(sp1,sp2),summarise,r=cor(ave.1,ave.2), sd1=sd(ave.1), sd2=sd(ave.2) ) 
sample.cnt.corr <- ddply(allpairs.cnt, .(sp1,sp2),summarise,r=cor(count.1,count.2) , sd1=sd(count.1), sd2=sd(count.2) ) 
sample.corr <- rbind(data.frame(var='ave',sample.ave.corr), data.frame(var='cnt',sample.cnt.corr))
x <- with(sample.corr, paste(sp1,sp2,sep='.'))
sample.corr$pair <- x

rescorr.end <- merge(res_corr, sample.corr)
rescorr.end$pair <- reorder(rescorr.end$pair, rescorr.end$r)
rescorr.end$var <- factor(rescorr.end$var, levels=c('ave','cnt'),labels =c('Mean', 'Total Count'))
rescorr.end$type <- factor(rescorr.end$type, levels=c('pairwise','tendim'), labels=c('Bivariate', 'Ten-Dimensional'))
colnames(rescorr.end)[4] <- 'Prior'
rescorr.end$Prior <- factor(rescorr.end$Prior, levels=c("iw","siw","ht","ss","iw.sc"),labels=c('IW','SIW', 'HIWht', 'BMMmu','IWsc'))

sample.var <- ddply(spcorr.dt, .(abbrev), summarise, ave = var(ave.add), cnt=var(count.add))
x <- with(sample.var, paste(abbrev,abbrev,sep='.'))
sample.var$pair <- x
sample.var <- melt(sample.var, id.vars=c('abbrev','pair'), variable_name='var')

resvar.end <- merge(res_var,sample.var)

resvar.end$var <- factor(resvar.end$var, levels=c('ave','cnt'),labels =c('Average', 'Total Count'))
resvar.end$type <- factor(resvar.end$type, levels=c('pairwise','tendim'), labels=c('Bivariate', 'Ten-Dimensional'))
colnames(resvar.end)[4] <- 'Prior'
resvar.end$Prior <- factor(resvar.end$Prior, levels=c("iw","siw","ht","ss","iw.sc"),labels=c('IW','SIW', 'HIWht', 'BMMmu','IWsc'))

resvar.end$mean2 <- resvar.end$mean
resvar.end$mean2[resvar.end$Prior=='IWsc'] <- resvar.end$mean2[resvar.end$Prior=='IWsc']*resvar.end$value[resvar.end$Prior=='IWsc']


# prior effect on correlations plots
pdf(file='report/figs/rescorr.pdf', height=4.5)
d <- subset(rescorr.end, Prior != 'IWsc')
qplot(data=d ,y=mean,x=r,shape=Prior, color=Prior,facets=var~type,ylab='Correlation posterior Mean', xlab='Pearson correlation coefficient') + geom_abline(1)+ theme(legend.position='bottom')
dev.off()

# prior effect on variance  plots
pdf(file='report/figs/resvar.pdf', height=4)
qplot(data=resvar.end ,y=mean2,x=value,shape=Prior, color=Prior, 
       ylab='Standard deviation posterior Mean', xlab='Sample standard deviation') + 
      facet_wrap(facets=~var,scales='free') + 
  geom_abline(1)+ theme(legend.position='bottom') + scale_x_log10()+ scale_y_log10()
dev.off()


# bird correaltions plots
pdf(file='report/figs/corrmat.pdf', height=5)
d1 <- subset(rescorr.end, Prior=='BMMmu' & type=='Ten-Dimensional' & var=='Average')
qplot(data=d1, sp1,sp2,geom='tile',fill=mean,xlab='',ylab='') +
scale_fill_gradient2(name=expression("Posterior Mean" * ~ rho),low='black',high='red',breaks=seq(-1, 1, by = 0.25), limits = c(-1, 1))+ 
theme(legend.position = 'top',panel.background = element_blank(), legend.direction = "horizontal")  +
  guides(fill = guide_colorbar(barwidth = 15, barheight = 1, title.position = "top", title.hjust = 0.5))
dev.off()
# Identify oven,revi and cswa are 3 correlated species,
#  lefl seems un-correlated with all the rest, specially with oven and amro

pdf(file='report/figs/spcor.pdf', height=5)
d <-  melt(all10.ave[, c('year','ave.CSWA','ave.REVI','ave.OVEN','ave.LEFL','ave.AMRO')],id.vars='year')
levels(d$variable) <- gsub('ave.', '', levels(d$variable))
colnames(d)[2] <- 'species'
qplot(data=d, x=year, y=value, color=species, geom='line', ylab='Difference respcet to historic mean') +
theme(legend.position='bottom')
dev.off()

# for discussion, inverse gamma plot
library(MCMCpack)
x <- seq(0.0001, 1/2,,1001)
ig <- dinvgamma(x, 1,1/2)
pdf(file='report/figs/ig.pdf', height=3.5)
qplot(x,ig, geom='line',ylab='density',xlab='') + geom_vline(xintercept=0.01, color=I('red')) + theme_bw()
dev.off()


#===============================================================
# 2) Run 6 'cells' of models, calling file codemodels.R to run them and save results. 
# mod.dd:     dif beta, dif sigma 
# mod.ds:     dif beta, same sigma
# mod.hs:     hier beta, same sigma
# res.dh.lin: dif beta, hier sigma, lin model, using jags 
# res.dh.q:   dif beta, hier sigma, quad model, using jags
# res.hd.lin: hier beta, dif sigma, lin model, using jags
# res.hd.q:   hier beta, dif sigma, quad model, using jags
# res.hh.lin: hier beta, hier sigma, lin model, using jags
# res.hh.q:   hier beta, hier sigma, quad model, using jags    

source('R\\codemodels.R') # jags text models are here

# runjags : function to run the jags model, and obtain the coda samples
# parameters: 
# d: data set, ave.add is the response
# model: the character object with the jags model
# l    : nbr of predictors, l=2 for linear and l=3 for quadratic
# lg   : if lg='logs' then the response is log transformed. 
runjags <- function(d, model, l,lg) {
  # d <- subset(bird.yeartotal, forest=='9020')
  if (lg=='logs') d$ave.add <- log(d$ave.add)
  m <- max(d$yearc)
  dend <- with(d, d[yearc==m, c('abbrev','ave.add')])
  dat = list(y = d$ave.add , 
             abbrev = as.numeric(d$abbrev) ,
             year= d$yearc, 
             n = nrow(d), 
             ns = nlevels(d$abbrev),
             R = diag(l),
             end = m+1,
             yend= dend$ave.add,
             abbvrev.end = dend$abbrev)
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=100)
  update(m, 1000)
  coda.samples(m, c('alpha','lambda','mu',"sigma.be",'beta','sigma.e','rate'), 3000)
}

# Create ml, a data.frame with the arguments needed to run each model. 
# this will be the argument of the runmod function. 
g  <- expand.grid('mtext' , c('dd', 'dh', 'ds', 'hd', 'hh', 'hs'), c('lin','quad'),  lg = c('logs', 'count') ) 
g$l <- 3 ; g$l[g$Var3=='lin'] <- 2
ml <- data.frame(mn = paste(g[,1],g[,2],g[,3],sep='.'), l=g$l, lg=g$lg)

runmod <- function( mdat ) { 
  m <- as.character(mdat$mn); ls <- mdat$l ; lgs <-as.character(mdat$lg)
  res  <-   dlply(bird.yeartotal, .(forestN), .fun=runjags, model=get(m), l=ls,lg=lgs)
}

# -----------------
# Run all 72 models
# results ia list of 24 model results, each element is a list of 
# 3 models one for each forest. 
ptm <- proc.time()
results <- alply(ml, .margins=1, runmod,  .progress='text')
proc.time() - ptm
#------------------

# to name results (I need to make it automatically) 
a <- adply(as.matrix(ml$mn),1, .fun=function(x) unlist(strsplit(x, '\\.')) )
a1 <- a[,3]
a2 <- paste(substr(a[,4], 1,1),substr(ml$lg, 1,1),sep='') 
names(results) <- paste('res',a1,a2, sep='.')

# Save the model results
save(results, file='modresults.Rdata')

#===============================================================
# 2 ) Summary of models, plots and tables 
load('modresults73.Rdata')

# function to collect: take the models 'mr' and get the 'par' parameters, col is the column in the summary
rescol <- function(mr, prm, col) {
  mds <- grep(mr , names(results))
  collection <- NULL
  for (i in mds) {
    d <- ldply(results[[i]], function(x) {
      cfn <- attributes(x[[1]])$dimnames[[2]]
      pr <- agrep(prm, cfn)
      s <- summary(x[,cfn[pr]])$quantile[,col]
      data.frame(estimate=s,row.names=NULL)
    })
    a <- unlist(strsplit(names(results)[i], '\\.'))[2]
    b <- unlist(strsplit(names(results)[i], '\\.'))[3]
    collection <-  rbind(collection,data.frame(mod=a,response=b, d))
  }
  return(collection)
}

# 2.1 Separate regresions: lm without random terms one per specie*forest, 
# response variable count, quadratic model 

sepcoef <- rescol('dd', 'beta 2', col=3)
sepsig.e <- rescol('dd', 'sigma.e', col=3)

qplot(data=subset(sepsig.e,response='ql'), y=estimate,x=forestN  )


pdf('figs/hist_m1.pdf')
#qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_wrap(facets=forestN~parameter,scale='free')
qplot(data=sepcoef, x=estimate,y=..density..,geom='histogram') +facet_wrap(facets=forestN~response,scale='free')
dev.off()

# bivariate plot
#x<- results$res.dd.ql$Superior
biv.dd <- ldply(results$res.dd.ql, function(x) {
  cfn <- attributes(x[[1]])$dimnames[[2]]
  dd <- data.frame( 
    b1 = summary( x[ ,cfn[agrep('beta 1',cfn)] ] )$quantile[,3],
    b2 = summary( x[ ,cfn[agrep('beta 2',cfn)] ] )$quantile[,3],
    b3 = summary( x[ ,cfn[agrep('beta 3',cfn)] ] )$quantile[,3],
    sig = summary( x[ ,cfn[grep('sigma.e',cfn)] ] )$quantile[,3])
  }
)
pdf(file='figs/scat_m1.pdf')
ggpairs(data=biv.dd, columns=2:5,color='forestN')
dev.off()

# beta joint dist  on hh.ql model
aux <- data.frame(rescol('hh.ql', prm=c('beta 2') ,col=3),
                  b3=rescol('hh.ql', prm=c('beta 3') ,col=3)$estimate)



# --------------------------

# 2.2 Diagnostic 
# gm compute the multivariate gelman.diag for one model  
gm <- function(x, multi=TRUE) {
  param <- attributes(x[[1]])$dimnames[[2]]
  aux.s <- grep('sigma.b', param)
  if (length(aux.s)==9) aux <- param[-aux.s[c(4,7,8)] ]
  if (length(aux.s)==4) aux <- param[-aux.s[3] ]
  if (length(aux.s)==0) aux <- param
  g <- gelman.diag(x[,aux])
  if (multi==TRUE)  out <-  g$mpsrf
  if (multi==FALSE) out <-  g
  return(out)
}

gr.multi <- ldply(results, function(z) ldply(z, gm) )
#subset(gr.multi, lg=='logs' & l==2)
#subset(gr.multi, V1 > 1.1)

gr.multi$dim <- ldply(results, function(z) ldply(z, function(x) length(attributes(x[[1]])$dimnames[[2]])) )$V1 

pdf('figs/gelman.pdf')  
qplot(data=gr.multi,x=1:72,ymin=1,ymax=V1,color=dim,geom='linerange', size=I(1),ylim=c(.95,1.3)) + geom_hline(yintercept=1.1, colour='red')
dev.off()


# individually check the models with big multivariate gelman
#mod <- runmod(ml[7,])
mod <- results[[7]]

cfn <- attributes(mod[[1]][[1]])$dimnames[[2]]
aux <- c(grep('alpha', cfn),grep('lambda', cfn) ) # grep('mu', cfn), agrep('sigma b', cfn)[c(1:3,5,6,9)])
llply(mod, function(x) gelman.diag(x[,grep('beta',cfn)]))
plot( mod[[1]][,c('alpha', 'lambda')], auto.layout=FALSE)
# compare simulations from alpha and lambda
as <- data.frame(as.matrix(mod[[2]][,c('alpha','lambda')]), 
                 chain=rep(1:3,each=3000))
qplot(data=as, x=lambda,y=alpha, size=I(.5),facets=chain~.)
qplot(data=as, x=rep(1:3000,3), y=lambda/alpha, size=I(.5),color=chain)
# rates ??? 
plot( mod[[1]][,c('rate[1]')], auto.layout=FALSE)
summary( mod[[1]][,grep('rate', cfn)])

#---------------------------------------
# 2.3 Slopes and rates plot to compare differents models

# collect slopes for linear models in logs
slopes.ll <- rescol('ll','beta 2')

# collect rates for quadratic models in logs 
rates.ql <- rescol('ql', 'rate') 

# extreme values in the rates: PROBLEM
ddply(rates.ql, .(mod, forestN), function(x) quantile(x$estimate, probs=seq(0,1,.1)))
subset(rates.ql, estimate > 200)

# function to plot 
# lms is a character having the two models to compare "a.b" 
# compare model a with model b. 
sloplot <- function(lms, coefs) {
  lms <- unlist(strsplit(lms, '\\.'))
  a <- lms[1] ; b <- lms[2]
  aux <- data.frame( subset(coefs, mod==a), subset(coefs, mod==b) )
  titu <- paste(a,b)
  qplot(data=aux, x=estimate, y=estimate.1, facets=~forestN, main=titu,xlab=a, ylab=b) + geom_abline(slope=1)
}
#sloplot("dd.dd")
plts <- list('dd.hd', 'dh.hh', 'ds.hs', 'hd.hh', 'hd.hs', 'hh.hs')
slp.plts <- llply(plts,  sloplot, coefs=slopes.ll)

aux <- subset(rates.ql, abs(estimate) <= 3)
rate.plts <- llply(plts,  sloplot, coefs=aux)

pdf('figs/slp_bet.pdf')
grid.arrange(slp.plts[[1]],slp.plts[[2]],slp.plts[[3]], nrow=3)
dev.off()

pdf('figs/rates.pdf')
grid.arrange(rate.plts[[1]],rate.plts[[2]],rate.plts[[3]], nrow=3)
dev.off()


#----------------------------------
# 2.4 Study posterior for mu and sigma.be
# we can study the sigma.be posterior in different context, 
# all hs, hh, and hd models has posterior for sigma.be


x <- results$res.hh.ql$Superior 

sig <- data.frame(as.matrix(x[,'mu[2]']),
                  as.matrix(x[,'mu[3]']),
                  as.matrix(x[,'sigma.be[2,3]']), 
                  as.matrix(x[,'sigma.be[3,3]']),
                  as.matrix(x[,'sigma.be[2,2]']) )
colnames(sig) <- c('mu2','mu3','s23', 's33', 's22')
rho23 <- adply(sig, .margins=1, summarise, rho=s23/sqrt(s33*s22))


pdf('quest.pdf')
qplot(data=aux, x=estimate,y=b3,color=forestN, main='Slope vs Quadratic on HH model', xlab='b2')
qplot(data=sig, x=mu2, y=mu3,size=I(.6),main='Posterior distribution for betas mean')+geom_density2d(size=I(1))
qplot(data=rho23, x=rho, geom='histogram', main='Posterior distribution for Correlation among b2 and b3')
dev.off()


covplot <- function(x) {
  param <- attributes(x[[1]])$dimnames[[2]]
  aux.s <- grep('sigma.b', param)
  l <- sqrt(length(aux.s))
  simx <- data.frame(as.matrix(x[[1]][,aux.s]))
  
  if (l==3) {
    c13 <- matrix(c(0,1,0,0,0,1,1,0,0), ncol=3)
    r13 <- matrix(c(0,0,1,1,0,0,0,1,0), ncol=3)
    mat <- list() 
    for (i in 1:dim(simx)[1])   {
      aux1 <- as.numeric(simx[i,])
      aux2 <- matrix(aux1,ncol=l,nrow=l)%*%c13
      aux <- r13%*%aux2
      mat[[i]] <- aux
    }
  }
  if (l==2) {
      mat <- list() 
        for (i in 1:dim(simx)[1])   {
              aux1 <- as.numeric(simx[i,])
                mat[[i]] <- matrix(aux,ncol=l,nrow=l)
              }
        }
  param = list(prob = 0.5, mat = mat)
  covdata =  VisCov("User defined distribution", param , title.logical = FALSE)  
  return(covdata)  
}

pdf(file='var_hhql.pdf')
cvd <- covplot(results$res.hd.ql$Chippewa)
dev.off()


postscript(file='var_hhll.ps')
covplot(results$res.hs.ql$Superior)
dev.off()


colnames(simx) <- c('s11', 's21','s12', 's22')
simx$rho <- with(simx,  s21/sqrt(s11*s22))
simx$sd11 <- log(sqrt(simx$s11))
simx$sd22 <- log(sqrt(simx$s22))
p1 <- qplot(data=simx, x=sd11, geom='histogram')
p2 <- qplot(data=simx, x=sd22, geom='histogram')
p3 <- qplot(data=simx, x=rho, geom='histogram')
p4 <- qplot(data=simx, x=sd22,y=sd11, size=I(.5))
p5 <- qplot(data=simx, x=rho,y=sd11, size=I(.5))
p6 <- qplot(data=simx, x=rho,y=sd22, size=I(.5))
grid.arrange(p1,p4,p2,p5,p3,p6, nrow=3)

plot(density(simx$sd22))

#--------------------------------------------
# Scaled IW

runjags <- function(d, model, l,lg) {
  # d <- subset(bird.yeartotal, forest=='9020')
  if (lg=='logs') d$ave.add <- log(d$ave.add)
  m <- max(d$yearc)
  dend <- with(d, d[yearc==m, c('abbrev','ave.add')])
  dat = list(y = d$ave.add , 
             abbrev = as.numeric(d$abbrev) ,
             year= d$yearc, 
             n = nrow(d), 
             ns = nlevels(d$abbrev),
             R = diag(l),
             end = m+1,
             yend= dend$ave.add,
             abbvrev.end = dend$abbrev)
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=100)
  update(m, 100)
  #coda.samples(m, c('alpha', ), 30)
  coda.samples(m, c('alpha','lambda','rho23','mu',"sigma.be",'beta','sigma.e'), 3000)
}


xdat <- subset(bird.yeartotal, forestN=='Superior')

ressup.SIW <- runjags(xdat, mtext.hh.siw, l=3, lg='logs')  
ressup.IW  <- runjags(xdat, mtext.hh.iw, l=3, lg='logs')  

param <- attributes(ressup.SIW[[1]])$dimnames[[2]]

betamed.siw <- data.frame(b1= summary(ressup.SIW[, param[agrep('beta 1', param)] ])$quantile[,3],
                          b2=summary(ressup.SIW[, param[agrep('beta 2', param)] ])$quantile[,3],
                          b3= summary(ressup.SIW[, param[agrep('beta 3', param)] ])$quantile[,3] )
qplot(data=betamed.siw, x=b2, y=b3)


# visualization for SIW model
simx <- data.frame(as.matrix(ressup.SIW[, c('sigma.be[2]','sigma.be[3]','rho23')]))

mat <- list() 
for (i in 1:dim(simx)[1])   {
  aux1 <- as.numeric(simx[i,])
  mat[[i]] <- matrix(aux1[c(1,3,3,2)],ncol=2)
}
dfp <- ldply(mat, function(x) (all(eigen(x)$values > 0)) )
summary(dfp) # there are 6 matrix not DF positive (check)

mat1 <- mat[!dfp$V1] 

param = list(prob = 0.5, mat = mat)
covdata =  VisCov("User defined distribution", param , title.logical = FALSE)  



# correlation b2 b3 using IW or SIW
save(ressup.IW, ressup.SIW, file='siw.Rdata')
summary(ressup.IW[, 'rho23'])
summary(ressup.SIW[, 'rho23'])

p <- c('sigma.be[1,1]','sigma.be[3,3]','sigma.be[2,2]','sigma.be[1,2]','sigma.be[1,3]','sigma.be[2,3]' )
summary(results$res.hh.ql$Superior[,p])

#========================================================
# shiny application to plot average over years 
#setwd('~\\GitHub\\mnfb_yearly\\data')
# runApp("shiny1")
  
# if (l==4) {
#   aux <- param[aux.s[c(1:2,4)]]
#   simx <- data.frame(as.matrix(x[[1]][,aux]))
#   mat <- list() 
#   for (i in 1:dim(simx)[1])   {
#     aux <- as.numeric(simx[i,c(1,2,2,3)])
#     mat[[i]] <- matrix(aux,ncol=2,nrow=2)
#   }
# } 
# if (length(aux.s)==9) {
#   aux <- param[aux.s[c(1:3, 5:6, 9)]]
#   simx <- data.frame(as.matrix(x[[1]][,aux]))
#   mat <- list() 
#   for (i in 1:dim(simx)[1])   {
#     aux <- as.numeric(simx[i,c(1,2,2,3)])
#     mat[[i]] <- matrix(aux,ncol=2,nrow=2)
#   }
# }
#===============

# Plot conditional distribution
rho_f = function(rho, v=3, d=2, sigma1=1, sigma2=1) {
  # identity matrix assumed for location parameter of IW distribution
  tr = (sigma1^2+sigma2^2)/(sigma1^2*sigma2^2*(1-rho^2))
  (1-rho^2)^(-(v+d+1)/2)*exp(-tr/2)
}
d = 2
int = integrate(function(x) rho_f(x,d+1,d),-1,1)
curve(rho_f(x,d+1,d)/int$value, -1, 1, col="red", n=10001)



# conditional distribution on sIW 

library(MCMCpack)

r <- numeric()
for (i in 1:20 ) {
  w1  <- rWishart(1, 3, diag(2))
  w  <- w1[ , ,1] 
  l  <- exp(rnorm(2, 0,1))
  r[i] <- prod(l)*w[1,2]
} 




