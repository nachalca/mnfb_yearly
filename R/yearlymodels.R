
# Modelling the yearly average of all species. 
# Data are created on mnfb repository, with the yearly_data.R code. 

bird.yeartotal <- read.csv('~\\GitHub\\mnfb\\data\\bird_yeartotal.csv')

# compute average using the count.add (this step should be on yearly_data.R)
bird.yeartotal$ave.add <- with(bird.yeartotal, count.add/samples)

setwd('~\\GitHub\\mnfb_yearly')

# libraries 
library(xtable)
library(plyr)
library(reshape2)
library(ggplot2)
library(lme4)
library(shiny)
library(GGally)


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

# 1) Simplest model: linear models without random terms one per specie*forest
m_f = function(x) {
  m = lm(log(ave.add) ~ I(year-2000) + I((year-2000)^2), x)
  data.frame(parameter=c("b0","b1","b2","sigma"),
             estimate=c(coef(m), summary(m)$sigma))
}

models.lm = ddply(bird.yeartotal, .(forestN,abbrev), m_f)

tab1 <- ddply(models.lm, .(forestN,parameter),function(x) summary(x$estimate) ) 
tab1 <- xtable(tab1, digits=3)
print(tab1, file="summ_modlm.tex", include.rownames=FALSE)

den.fun <- function(z) {
  d <- density(z$estimate)
  data.frame(knots=d$x,dens=d$y) 
}
densities <- ddply(models.lm, .(forestN, parameter), den.fun)

postscript('figs/hist_m1.ps')
#qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_grid(facets=parameter~forest,scale='free')
qplot(data=densities, x=knots,y=dens,geom='line')+facet_wrap(facets=parameter~forestN, scales="free",ncol=3)  
dev.off()

aux <- reshape(models.lm , direction='wide', timevar='parameter', idvar=c('forestN','abbrev') )
postscript('figs/scat_m1.ps')
ggpairs(data=aux, columns=3:6, color='forestN')
dev.off()

#=====================================================
# 2) Include random terms in the model

# model using lmer, adding random effects. 

mrnd_f <- function(x) {
  x$year <- x$year-1994
  x$abbrev <- factor(x$abbrev)
  x$lave <- log(x$ave) 
  x$lave[x$lave==-Inf] <- NA
  
  mod <-lmer(lave ~ year+I(year^2)+(year+I(year^2)|abbrev),data=x)  
  
  data.frame(parameter=c("b0","b1","b2","residual"),
            mean=c(fixef(mod)[1:3],0), variance=c(attributes(VarCorr(mod)$abbrev)$stddev,attributes(VarCorr(mod))$sc) )
}

models.rnd = ddply(bird.yeartotal, .(forest), mrnd_f)
tab2 <- xtable(models.rnd, digits=3)
print(tab2, file="summ_modrnd.tex", include.rownames=FALSE)

mrnd_f2 <- function(x) {
  x$year <- x$year-1994
  x$abbrev <- factor(x$abbrev)
  x$lave <- log(x$ave) 
  x$lave[x$lave==-Inf] <- NA
  
  mod <-lmer(lave ~ year+I(year^2)+(year+I(year^2)|abbrev),data=x)  
  data.frame(coef(mod)$abbrev)
}
models.rndeff = ddply(bird.yeartotal, .(forest), mrnd_f2)
aux <- melt(data=models.rndeff, id.vars='forest')

postscript('..\\figs\\box_m2.pdf')
qplot(data=aux, x=factor(forest), group=factor(forest), y=value, geom='boxplot')+facet_grid(facets=variable~., scale='free')
dev.off()

# -----------------------------------------------------
# -----------------------------------------------------

#3) shiny application to plot average over years 
runApp("shiny1")

