
# Modelling the yearly average of all species. 
# Data are created on mnfb repository, with the yearly_data.R code. 

bird.yeartotal <- read.csv('~\\GitHub\\mnfb\\data\\bird_yeartotal.csv')

# compute average using the count.add (this step should be on yearly_data.R)
bird.yeartotal$ave.add <- with(bird.yeartotal, count.add/samples)

# create a centered year variable
bird.yeartotal$yearc <- bird.yeartotal$year-2000
setwd('~\\GitHub\\mnfb_yearly')

# libraries 
library(xtable)
library(plyr)
library(reshape2)
library(lme4)
library(shiny)
library(ggplot2)
library(GGally)
library(rjags)
library(gridExtra)

#-------------------------------------------------
# 0) data description
tab_f <- function(x,sp) {
  x <- x[x$abbrev %in% sp,]
  data.frame(x[,c('count','abbrev')])
}

# get the 15 more abundant species
aux <- sort(with(bird.yeartotal, tapply(count,abbrev,sum)),decreasing=T )[1:15]
spn <- names(aux)
tab <- ddply(bird.yeartotal, .(year,forestN), tab_f, sp=spn)

# pull out year 2007 to compare with the online report
aux <- reshape(subset(tab, year==2007) , direction='wide', timevar='forestN', idvar='abbrev',drop='year')
aux <- aux[order(aux[,2],decreasing=T),]
colnames(aux) <- c('Specie',levels(bird.yeartotal$forestN) )
tab = xtable(aux)
print(tab, file="highcounts.tex", include.rownames=FALSE)

# Overall trend within forest 
totales <- ddply(bird.yeartotal, .(year,forestN), summarise, count=sum(count) )
totales$forest <- factor(totales$forest)

postscript('figs/rawtrend.ps')
qplot(data=totales,x=year,xlab='Year',ylab='Bird Count',y=count,shape=forestN,color=forestN, size=I(3))+geom_line()
dev.off()
#===============================================================

# 1) Run 6 'cells' of models, calling file codemodels.R to run them and save results. 
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
  dat = list(y = d$ave.add , 
             abbrev = as.numeric(d$abbrev) ,
             year= d$yearc, 
             n = nrow(d), 
             ns = nlevels(d$abbrev),
             R = diag(l))
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=100)
  update(m, 1000)
  coda.samples(m, c('alpha','lambda','mu',"sigma.be",'beta','sigma.e'), 3000)
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
load('modresults.Rdata')

# 2.1 Separate regresions: lm without random terms one per specie*forest, 
# response variable count, quadratic model 
pr <- agrep( "beta 2", pnam.dd.cc)
mr <- grep('dd', names(results))
sepcoef <- NULL
for (i in mr) {
  d <- ldply(results[[i]], function(x) {
    s <- summary(x[,pnam.dd.cc[pr]])$quantile[,3]
  data.frame(estimate=s,row.names=NULL)
  })
  sepcoef <-  rbind(sepcoef,data.frame(mod=names(results)[i], d))
}
postscript('figs/hist_m1.ps')
#qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_wrap(facets=forestN~parameter,scale='free')
qplot(data=sepcoef, x=estimate,y=..density..,geom='histogram') +facet_wrap(facets=forestN~mod,scale='free')
dev.off()

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

subset(gr.multi, lg=='logs' & l==2)
subset(gr.multi, V1 > 1.1)

# individually check the models with big multivariate gelman
# mod <- runmod(ml[11,])
mod <- results[[20]]
cfn <- attributes(mod[[1]][[1]])$dimnames[[2]]
aux <- c(grep('alpha', cfn),grep('lambda', cfn) ) # grep('mu', cfn), agrep('sigma b', cfn)[c(1:3,5,6,9)])
llply(mod, function(x) gelman.diag(x[,aux]))

plot( mod[[2]][,c('alpha', 'lambda')], auto.layout=FALSE)

# compare simulations from alpha and lambda
as <- data.frame(as.matrix(mod[[2]][,c('alpha','lambda')]), 
                 chain=rep(1:3,each=3000))
qplot(data=as, x=lambda,y=alpha, size=I(.5),facets=chain~.)
qplot(data=as, x=rep(1:3000,3), y=lambda/alpha, size=I(.5),facets=chain~.)


#---------------------------------------
# 2.3 Slopes plot to compare differents models

# collect slopes for linear models in logs
mr <- grep('ll', names(results))
slopes.ll <- NULL
for (i in mr) {
  d <- ldply(results[[i]], function(x) {
    cfn <- attributes(x[[1]])$dimnames[[2]]
    pr <- agrep('beta 2', cfn)
    s <- summary(x[,cfn[pr]])$quantile[,3]
    data.frame(estimate=s,row.names=NULL)
  })
  a <- unlist(strsplit(names(results)[i], '\\.'))[2]
  slopes.ll <-  rbind(slopes.ll,data.frame(mod=a, d))
}

# function to plot 
# lms is a character having the two models to compare "a.b" 
# compare model a with model b. 
sloplot <- function(lms) {
  lms <- unlist(strsplit(lms, '\\.'))
  a <- lms[1] ; b <- lms[2]
  aux <- data.frame( subset(slopes.ll, mod==a), subset(slopes.ll, mod==b) )
  titu <- paste(a,b)
  qplot(data=aux, x=estimate, y=estimate.1, facets=~forestN, main=titu,xlab=a, ylab=b) + geom_abline(slope=1)
}
#sloplot("dd.dd")
plts <- list('dd.hd', 'dh.hh', 'ds.hs', 'hd.hh', 'hd.hs', 'hh.hs')
slp.plts <- llply(plts,  sloplot)

postscript('figs/slp_bet.ps')
grid.arrange(slp.plts[[1]],slp.plts[[2]],slp.plts[[3]], nrow=3)
dev.off()



#========================================================
# shiny application to plot average over years 
runApp("shiny1")

