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
