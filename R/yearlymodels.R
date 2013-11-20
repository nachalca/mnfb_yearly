
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
#--------------

# 2) Run 6 'cells' of models, calling file codemodels.R to run them and save results. 

# mod.dd:     dif beta, dif sigma, using lm()
# mod.ds:     dif beta, same sigma, using lm() 
# mod.hs:     hier beta, same sigma, using lmer()
# res.dh.lin: dif beta, hier sigma, lin model, using jags 
# res.dh.q:   dif beta, hier sigma, quad model, using jags
# res.hd.lin: hier beta, dif sigma, lin model, using jags
# res.hd.q:   hier beta, dif sigma, quad model, using jags
# res.hh.lin: hier beta, hier sigma, lin model, using jags
# res.hh.q:   hier beta, hier sigma, quad model, using jags    

source('R\\codemodels.R')
save(mod.dd, mod.ds, mod.hs, res.dh.lin, 
     res.dh.q, res.hd.lin,res.hd.q, res.hh.lin,res.hh.q, 
     file='modresults.Rdata')
#===============================================================
# 3 ) Summary of models, plots and tables 

# Separate regresions: lm without random terms one per specie*forest
qplot(data=subset(mod.dd, parameter=='b1'), x=model, y=aic,color=ci.sig, facets=~forestN)
with(mod.dd, table(parameter, ci.sig,model))
tab1 <- ddply(mod.dd, .(forestN,model,parameter),function(x) quantile(x$estimate, probs=c(0.25, 0.5, 0.975)) ) 
tab1 <- xtable(tab1, digits=3, caption='Percentiles for Separate Regresion model coefficients')
print(tab1, file="summ_modlm.tex", include.rownames=FALSE)

postscript('figs/hist_m1.ps')
#qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_wrap(facets=forestN~parameter,scale='free')
qplot(data=subset(mod.dd, parameter=='b1'), x=estimate,y=..density..,geom='histogram') +facet_wrap(facets=forestN~model,scale='free')
dev.off()

aux <- reshape(mod.dd , direction='wide', timevar='parameter', idvar=c('forestN','abbrev') )
postscript('figs/scat_m1.ps')
ggpairs(data=aux, columns=3:6, color='forestN')
dev.off()

levels(models.rnd$parameter)[c(1:2,4,5)] <- c('b0','sigma', 'b1', 'b2')
tab2 <- xtable(models.rnd, digits=3, caption='Random coefficent regresion results')
print(tab2, file="summ_modrnd.tex", include.rownames=FALSE)


# Bayesian models diag, models 4,5,6

par.names <- attributes(res.wishart$Chippewa[[1]])$dimnames[[2]]
par.names.lin <- attributes(res.wishart.lin$Chippewa[[1]])$dimnames[[2]]
# gelman diagnostic, only 6 param from the covariance matrix  
gquad <- llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names[c(1,221:227,229:230,233)] ]) )
llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names.lin[c(1,148:152,154)] ]) )

#plot(res.wishart[[1]][, par.names[1:3]])
# rearrange simulations for ploting
resdd <- function(res) {
  p  <- dim(res[[1]])[2]
  it <- dim(res[[1]])[1]
  cn <- length(res)
  df <- NULL 
  for (i in 1:cn) {
    df1 <- data.frame(iter=1:it, chain=i, as.data.frame(res[[i]]))
    df <- rbind(df, df1)
  }
  melt(df, id.vars=c('iter', 'chain'), variable_name='par')    
}
dd <- ldply(res.wishart, resdd)
levels(dd$par)
# trace and den for alpha, lambda
dd1 <- subset(dd, par %in% c('alpha','lambda') )
qplot(data=dd1, x=iter, y=value,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
plot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')
# trace and den for mus
dd1 <- subset(dd, par %in% c('mu.1.','mu.2.','mu.3.') )
qplot(data=dd1, x=iter, y=value,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
qplot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')
# trace and den for vars
dd1 <- subset(dd, par %in% c('sigma.be.1.1.','sigma.be.2.2.','sigma.be.3.3.') )
qplot(data=dd1, x=iter, y=value,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
qplot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')
# trace and den for COV
dd1 <- subset(dd, par %in% c('sigma.be.1.2.','sigma.be.1.3.','sigma.be.2.3.') )
qplot(data=dd1, x=iter, y=value,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
qplot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')
tab3 <- ddply(dd, .(forestN, par), function(x) quantile(x$value, probs=c(.025, .5, .975)) )

t <- tab3[tab3$par %in% c('mu.1.','mu.2.','mu.3.','sigma.be.1.1.','sigma.be.2.2.','sigma.be.3.3.', 'sigma.be.2.1.','sigma.be.3.1.','sigma.be.2.3.'), ]
t <- xtable(t, digits=3, caption='Mean and variances for Wishar Prior model')
print(t, file="summ_wishart.tex", include.rownames=FALSE)
#========================================================


# collect results for bayesian model
summfun <- function(res) {
  s <- summary(res)
  data.frame(par = rownames(s$quantiles), s$quantiles[,c(1,3,5)] )
}

sum.wis <- ldply(res.wishart, summfun )
sum.wis.lin <- ldply(res.wishart.lin, summfun )


# With the results from all 'only linear' models compare slopes 

wishcoef <- ddply(sum.wis.lin, .(forestN), function(x) x[seq(3,147,2),-c(1,2)] )

sepcoef    <- subset(models.lm, parameter=='b1' & model=='llin')[,c(1,5,8)]
sepcoef$lo <- with(sepcoef, estimate - 1.96*sqrt(var) )                  
sepcoef$up <- with(sepcoef, estimate + 1.96*sqrt(var) )

rndcoef <- ddply(subset(models.rnd, .id=='lin'), .(forestN), function(x) x[-c(1:77),-c(1,2)] )
rndcoef$lo2 <- with(rndcoef, mean - 1.96*sqrt(variance) )                  
rndcoef$up2 <- with(rndcoef, mean + 1.96*sqrt(variance) )

# Plot
aux <- data.frame(sepcoef, rndcoef)
p1 <- qplot(data=aux, x=estimate, y=mean, facets=~forestN, main='Separate Regresions vs Random coef', xlab='sepreg', ylab='rndcoef') + geom_abline(slope=1)
p1e <- p1 + geom_errorbar(aes(ymax =up2, ymin=lo2), width = 0.05) 
p1e <- p1e + geom_errorbarh(aes(xmax=up, xmin=lo), width = 0.05) 


aux <- data.frame(sepcoef, wishcoef)
p2 <- qplot(data=aux, x=estimate, y=X50., facets=~forestN, main='Wishart Prior vs Separate Regresions', xlab='sepreg', ylab='wishart')  + geom_abline(slope=1)
p2e <- p2 + geom_errorbar(aes(ymax =X97.5., ymin=X2.5.), width = 0.05) 
p2e <- p2e + geom_errorbarh(aes(xmax=up, xmin=lo), width = 0.05) 

aux <- data.frame(rndcoef, wishcoef)
p3 <- qplot(data=aux, x=mean, y=X50., facets=~forestN, main='Random coef vs Wishart Prior', xlab='rndcoef', ylab='wishart')  + geom_abline(slope=1)
p3e <- p3 + geom_errorbar(aes(ymax =X97.5., ymin=X2.5.), width = 0.05) 
p3e <- p3e + geom_errorbarh(aes(xmax=up2, xmin=lo2), width = 0.05) 

grid.arrange(p1, p2, p3, nrow=3)

grid.arrange(p1e, p2e, p3e, nrow=3)

# or..... 
aux1 <- data.frame(sepcoef, rndcoef, compare='Random vs Separate')
aux2 <- data.frame(sepcoef, wishcoef, compare='Separate vs Wishart')
aux3 <- data.frame(rndcoef, wishcoef, compare='Random vs Wishart')

colnames(aux1)[8]   <- 'estimate2'
colnames(aux2)[7:9] <- c('lo2', 'estimate2', 'up2')
colnames(aux3)[c(3,5,6,8:10)] <- c('estimate','lo', 'up', 'lo2', 'estimate2', 'up2')

aux <- rbind(aux1[,c(1:2,4,5,8,10:12)],
      aux2[,c(1,2,4,5,8,7,9,10)],
      aux3[,c(1,3,5,6,9,8,10,11)]
      )
p <- qplot(data=aux, x=estimate, y=estimate2, facets=compare~forestN)
pe <- p + geom_errorbar(aes(ymax =up2, ymin=lo2), width = 0.05) 
pe <- pe + geom_errorbarh(aes(xmax=up, xmin=lo), width = 0.05) 

# or .... 
aux$est.sig <- with(aux, lo*up > 0)
aux$est2.sig <- with(aux, lo2*up2 > 0)


qplot(data=aux, x=estimate, y=estimate2,color=est.sig, shape=est2.sig, facets=compare~forestN)



#========================================================
# shiny application to plot average over years 
runApp("shiny1")

