
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

# 1) Separate regresions: lm without random terms one per specie*forest

m_f = function(x) {
  m <- list(4)
  m[[1]] = lm(ave.add ~ I(yearc) + I((yearc)^2), x) # quad 
  m[[2]] = lm(log(ave.add) ~ I(yearc) + I((yearc)^2), x) # lquad
  m[[3]] = lm(ave.add ~ I(yearc), x) # lin
  m[[4]] = lm(log(ave.add) ~ I(yearc), x) # llin
  names(m)  <- c('quad', 'lquad', 'lin', 'llin')
  
  ldply(m, function(ms) {
    data.frame(parameter=c(attributes(coef(ms))$names,"sigma"),
             estimate=c(coef(ms), summary(ms)$sigma),
             ci.sig = c( apply(confint(ms),1,prod) > 0 ,NA),aic=AIC(ms))
  }
  )
}
models.lm = ddply(bird.yeartotal, .(forestN,abbrev), m_f)
colnames(models.lm)[3] <- 'model'
levels(models.lm$parameter) <- c('b0', 'b2', 'b1', 'sigma')


qplot(data=subset(models.lm, parameter=='b1'), x=model, y=aic,color=ci.sig, facets=~forestN)

with(models.lm, table(parameter, ci.sig,model))


tab1 <- ddply(models.lm, .(forestN,model,parameter),function(x) quantile(x$estimate, probs=c(0.25, 0.5, 0.975)) ) 
tab1 <- xtable(tab1, digits=3, caption='Percentiles for Separate Regresion model coefficients')
print(tab1, file="summ_modlm.tex", include.rownames=FALSE)

postscript('figs/hist_m1.ps')
#qplot(data=models.lm, x=estimate,y=..density..,geom='histogram')+facet_wrap(facets=forestN~parameter,scale='free')
qplot(data=subset(models.lm, parameter=='b1'), x=estimate,y=..density..,geom='histogram') +facet_wrap(facets=forestN~model,scale='free')

dev.off()

aux <- reshape(models.lm , direction='wide', timevar='parameter', idvar=c('forestN','abbrev') )
postscript('figs/scat_m1.ps')
ggpairs(data=aux, columns=3:6, color='forestN')
dev.off()

#=====================================================
# 2) Include random terms in the model

# model using lmer, adding random effects. 

mrnd_f <- function(x) {
  x$lave <- log(x$ave.add) 
  m <- list(2)
  m[[1]] <-lmer(lave ~ yearc +(yearc|abbrev),data=x)
  m[[2]] <-lmer(lave ~ yearc+I(yearc^2)+(yearc+I(yearc^2)|abbrev),data=x)  
  names(m)  <- c('lin', 'quad')
  
  ldply(m, function(ms) {
    aux1 <- triu( attributes(VarCorr(ms)$abbrev)$correlation,1 ) 
    cor.val <- as.numeric(aux1[aux1!=0])  
    cor.name <- c('rho_01','rho_02','rho_12')[1:length(cor.val)]
    data.frame(parameter=c(names(fixef(ms)),"residual",cor.name),
             mean=c(fixef(ms),0, cor.val),
             variance=c(attributes(VarCorr(ms)$abbrev)$stddev,
                        attributes(VarCorr(ms))$sc,rep(0,length(cor.val))) )
  }
  )
}

#x <- subset(bird.yeartotal, forest==9020) 
#qplot(data=x, year,lave)
#qplot(data=x[x$abbrev=='ALFL',], year,lave) + geom_line()

models.rnd = ddply(bird.yeartotal, .(forestN), mrnd_f)

ci_fun <- function(x,al) {
    x <- as.numeric(x)
    ci <- qnorm(c(al/2, 1-al/2), mean=x[1], sd=x[2]^2)
    prod(ci) > 0
}
models.rnd$ci.sig<-  apply(models.rnd[, 4:5], 1, ci_fun, al=.05)
levels(models.rnd$parameter)[c(1:2,4,5)] <- c('b0','sigma', 'b1', 'b2')

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
  beta[1:3,j]   ~ dmnorm(mu,prec.be )
  eta.e[j]      ~ dgamma(alpha,lambda)
  sigma.e[j] <- 1/sqrt(eta.e[j])      
}

# hyperpriors
prec.be   ~ dwish(R, df)
sigma.be  <- inverse(prec.be)
for (i in 1:3) { mu[i] ~ dnorm(0,0.001) }

df     <- 4  
alpha  ~ dunif(0,100)
lambda ~ dunif(0,100)

# predictives 
beta.pred ~ dmnorm(mu,prec.be )
eta.epred ~ dgamma(alpha,lambda)
sigma.epred <- 1/sqrt(eta.epred)
}
"
model.wishart.lin = "
model {
for (i in 1:n) { 
y[i] ~ dnorm(eff[i]+slop[i] , eta.e[abbrev[i]]  )
eff[i]  <-  beta[1, abbrev[i]]
slop[i] <-  beta[2, abbrev[i]]*year[i]
}

# Priors.  
for (j in 1:ns) {
beta[1:2,j]   ~ dmnorm(mu,prec.be )
eta.e[j]      ~ dgamma(alpha,lambda)
sigma.e[j] <- 1/sqrt(eta.e[j])      
}

# hyperpriors
prec.be   ~ dwish(R, df)
sigma.be  <- inverse(prec.be)
for (i in 1:2) { mu[i] ~ dnorm(0,0.001) }

df     <- 4  
alpha  ~ dunif(0,100)
lambda ~ dunif(0,100)

# predictives 
beta.pred ~ dmnorm(mu,prec.be )
eta.epred ~ dgamma(alpha,lambda)
sigma.epred <- 1/sqrt(eta.epred)
}
"

runjags <- function(d, model, l) {
  # d <- subset(bird.yeartotal, forest=='9020')
  dat = list(y = log(d$ave.add)  , 
           abbrev = as.numeric(d$abbrev) ,
           year= d$yearc, 
           n = nrow(d), 
           ns = nlevels(d$abbrev),
           R = diag(l))
  m = jags.model(textConnection(model), dat, n.chains=3, n.adapt=1000)
  update(m, 1000)
  coda.samples(m, c('alpha','lambda','mu',"sigma.be",'beta.pred','sigma.epred'), 4000)
}


res.wishart <- dlply(bird.yeartotal, .(forestN), .fun=runjags, model=model.wishart, l=3)
res.wishart.lin <- dlply(bird.yeartotal, .(forestN), .fun=runjags, model=model.wishart.lin, l=2)

par.names <- attributes(res.wishart$Chippewa[[1]])$dimnames[[2]]

par.names.lin <- attributes(res.wishart.lin$Chippewa[[1]])$dimnames[[2]]

# gelman diagnostic, only 6 param from the covariance matrix  
llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names[c(1:11,13,14,17,18)] ]) )

llply(res.wishart, .fun=function(x) gelman.diag(x[, par.names.lin[-9] ]) )


plot(res.wishart[[1]][, par.names[1:3]])

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
qplot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')

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

# predictives for betas and sigma.e
dd1 <- subset(dd, par %in% c("beta.pred.1.","beta.pred.2.","beta.pred.3.","sigma.epred") )
qplot(data=dd1, x=iter, y=value,size=I(1), color=factor(chain))+ facet_grid(facets=par~forestN,scale='free_y')
qplot(data=dd1, x=value, y=..density..,geom='histogram')+ facet_wrap(facets=forestN~par,scale='free')


tab3 <- ddply(dd, .(forestN, par), function(x) quantile(x$value, probs=c(.025, .5, .975)) )
t <- tab3[tab3$par %in% c('mu.1.','mu.2.','mu.3.','sigma.be.1.1.','sigma.be.2.2.','sigma.be.3.3.', 'sigma.be.2.1.','sigma.be.3.1.','sigma.be.2.3.'), ]
t <- xtable(t, digits=3, caption='Mean and variances for Wishar Prior model')
print(t, file="summ_wishart.tex", include.rownames=FALSE)








# shiny application to plot average over years 
runApp("shiny1")

