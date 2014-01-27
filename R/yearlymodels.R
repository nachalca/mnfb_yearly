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
library(VisCov)
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
aux <- aux[order(aux[,2],decreasing=T),]
colnames(aux) <- c('Specie',levels(bird.yeartotal$forestN) )
tab = xtable(aux, caption='Total counts on 2007 for the 10 more abundant species')
print(tab, file="highcounts.tex", include.rownames=FALSE, caption.placement='top')

# Overall trend within forest 
totales <- ddply(bird.yeartotal, .(year,forestN), summarise, count=sum(count) )
totales$forest <- factor(totales$forest)

pdf('figs/rawtrend.pdf')
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


