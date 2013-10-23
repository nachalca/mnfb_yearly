
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

# 1) Simplest model: linear models without random terms one per specie*forest
m_f = function(x) {
  m = lm(log(ave.add) ~ I(yearc) + I((yearc)^2), x)
  data.frame(parameter=c("b0","b1","b2","sigma"),
             estimate=c(coef(m), summary(m)$sigma))
}

models.lm = ddply(bird.yeartotal, .(forestN,abbrev), m_f)

tab1 <- ddply(models.lm, .(forestN,parameter),function(x) quantile(x$estimate, probs=c(0.25, 0.5, 0.975)) ) 
tab1 <- xtable(tab1, digits=3, caption='Percentiles for Separate Regresion model coefficients')

print(tab1, file="summ_modlm.tex", include.rownames=FALSE)

#den.fun <- function(z) {
#  d <- density(z$estimate)
#  data.frame(knots=d$x,dens=d$y) 
#}
#densities <- ddply(models.lm, .(forestN, parameter), den.fun)

postscript('figs/hist_m1.ps')
qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_wrap(facets=forestN~parameter,scale='free')
#qplot(data=densities, x=knots,y=dens,geom='line')+facet_wrap(facets=parameter~forestN, scales="free_y",ncol=3)  
dev.off()

aux <- reshape(models.lm , direction='wide', timevar='parameter', idvar=c('forestN','abbrev') )
postscript('figs/scat_m1.ps')
ggpairs(data=aux, columns=3:6, color='forestN')
dev.off()

#=====================================================
# 2) Include random terms in the model

# model using lmer, adding random effects. 

mrnd_f <- function(x) {
  x$abbrev <- factor(x$abbrev)
  x$lave <- log(x$ave.add) 
  mod <-lmer(lave ~ yearc+I(yearc^2)+(yearc+I(yearc^2)|abbrev),data=x)  
  data.frame(parameter=c("b0","b1","b2","residual",'rho_01','rho_02','rho_12'),
             mean=c(fixef(mod)[1:3],0,attributes(VarCorr(mod)$abbrev)$correlation[c(2,3,6)]),
             variance=c(attributes(VarCorr(mod)$abbrev)$stddev,
                        attributes(VarCorr(mod))$sc,rep(0,3)) )
}

#x <- subset(bird.yeartotal, forest==9020) 
#qplot(data=x, year,lave)
#qplot(data=x[x$abbrev=='ALFL',], year,lave) + geom_line()

models.rnd = ddply(bird.yeartotal, .(forestN), mrnd_f)
tab2 <- xtable(models.rnd, digits=3, caption='Random coefficent regresion results')
print(tab2, file="summ_modrnd.tex", include.rownames=FALSE)

#=====================================================
# 3) Bayes model

# bayesian models: one model for each forest
# using wishart distribution for coviariance matrix of beta coef 

# n is total of years within each species (17) 
# ns is the total of species (73)

model.wishart = "
model {
for (i in 1:n) { 
y[i] ~ dnorm(eff[i]+slop[i]+quad[i] , eta.e[abbrev[i]]  )
eff[i]  <-  beta[1, abbrev[i]]
slop[i] <-  beta[2, abbrev[i]]*year[i]
quad[i] <-  beta[3, abbrev[i]]*year[i]^2
}

# Priors.  
for (j in 1:ns) {
  beta[1:3,j]   ~ dmnorm(mu,Sigbe )
  eta.e[j]       <- 1/sigma.e[j]^2
  sigma.e[j]      ~ dgamma(alpha,lambda)
}

Sigbe ~ dwish(R, df)
#sigma.be  <- inverse(Asigma.be[,])

for (i in 1:3) { mu[i] ~ dnorm(0,0.001) }

df     <- 4  
alpha  ~ dunif(0,100)
lambda ~ dunif(0,100)
}
"

runjags <- function(d, model) {
  # d <- subset(bird.yeartotal, forest=='9020')
  dat = list(y = log(d$ave.add)  , 
           abbrev = as.numeric(d$abbrev) ,
           year= d$yearc, 
           n = nrow(d), 
           ns = length(levels(d$abbrev)),
           R = diag(3))
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=500)
  coda.samples(m, c('alpha','lambda','mu',"Sigbe"), 4000)
}

res.wishart <- dlply(bird.yeartotal, .(forestN), .fun=runjags, model=model.wishart)

#quema <- function(x,st,end)   
#res.wishart <- llply(res.wishart1, quema, st=2001,end=4000 )
#aux <- window.mcmc(res.wishart$Chippewa, start=10, end=100)

par.names <- attributes(res.wishart$Chippewa[[1]])$dimnames[[2]]

# gelman diagnostic 
llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names[10:14] ]) )
# error ? when I included all 1 to 9 covariance parameters 
llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names[7:9] ]) )

resdd <- function(res,name,cn,it) {
  p <- length(name)
data.frame(
iter =  rep(1:it, times=p),
sim  =  unlist(res[, name ] ),
chain = rep(1:cn, each=p*it),
par  =  rep(rep(name,each=it),times=cn)
)
}

dd <- ldply(res.wishart, resdd, name=par.names[c(12:14,1,5,9,2,3,8) ], cn=3, it=1000)

postscript('figs/trace_wis.ps')
qplot(data=dd, x=iter, y=sim,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
dev.off()

tab3 <- ddply(dd, .(forestN, par), function(x) quantile(x$sim, probs=c(.025, .5, .975)) )
tab3 <- xtable(tab3, digits=3, caption='Wishar Prior model')
print(tab3, file="summ_wishart.tex", include.rownames=FALSE)


qplot(data=dd, x=sim, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')














# shiny application to plot average over years 
runApp("shiny1")

