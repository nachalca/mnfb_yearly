

# jags code for models models 

# inverse wishart
mtext.hh.iw <-  "
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

# use scale chi param
eta.e[j]   ~ dgamma(alpha/2, alpha*lambda/2)
sigma.e[j] <- 1/sqrt(eta.e[j])
}

# hyperpriors
prec.be   ~ dwish(R, df)
sigma.be  <- inverse(prec.be)
for (i in 1:3) { mu[i] ~ dnorm(0,0.001) }

df     <- 4  
alpha  ~ dunif(0, 1000)
lambda ~ dunif(0, 1000)
rho23 <- sigma.be[3,2]/sqrt(sigma.be[3,3]*sigma.be[2,2])

#predictives 
for (i in 1:ns) { 
ynext[i] ~ dnorm(e[i]+s[i]+q[i] , eta.e[abbvrev.end[i]])
e[i]  <-  beta[1, abbvrev.end[i]]
s[i] <-  beta[2, abbvrev.end[i]]*end
q[i] <-  beta[3, abbvrev.end[i]]*end^2
rate[i] <- ynext[i]/yend[abbvrev.end[i]] - 1 
}
}
"

# Scaled Inverse Wishart prior
mtext.hh.siw = "
model {
for (i in 1:n) { 
y[i] ~ dnorm(eff[i]+slop[i]+quad[i] , eta.e[abbrev[i]]  )
eff[i]  <-  beta[1, abbrev[i]]
slop[i] <-  beta[2, abbrev[i]]*year[i]
quad[i] <-  beta[3, abbrev[i]]*year[i]^2
}

# Priors.  
for (j in 1:ns) {
beta[1,j] <- xi[1]*beta.raw[1,j]
beta[2,j] <- xi[2]*beta.raw[2,j]
beta[3,j] <- xi[3]*beta.raw[3,j]
beta.raw[1:3,j]   ~ dmnorm(mu.raw, tau.raw)

# use scale chi param
eta.e[j]   ~ dgamma(alpha/2, alpha*lambda/2)
sigma.e[j] <- 1/sqrt(eta.e[j])
}

# hyperpriors
tau.raw   ~ dwish(R, df)
sigma.raw  <- inverse(tau.raw)

for (i in 1:3) { 
mu.raw[i] ~ dnorm(0,0.001) 
xi[i] ~ dunif(0, 100)
mu[i] <- xi[i]*mu.raw[i]
sigma.be[i] <- xi[i]*sqrt(sigma.raw[i,i])
}
rho12 <- xi[1]*xi[2]*sigma.raw[1,2]/(sigma.be[1]*sigma.be[2])
rho13 <- xi[1]*xi[3]*sigma.raw[1,3]/(sigma.be[1]*sigma.be[3])
rho23 <- xi[3]*xi[2]*sigma.raw[3,2]/(sigma.be[3]*sigma.be[2])

df     <- 4  
alpha  ~ dunif(0, 1000)
lambda ~ dunif(0, 1000)

#predictives 
#for (i in 1:ns) { 
#ynext[i] ~ dnorm(e[i]+s[i]+q[i] , eta.e[abbvrev.end[i]])
#e[i]  <-  beta[1, abbvrev.end[i]]
#s[i] <-  beta[2, abbvrev.end[i]]*end
#q[i] <-  beta[3, abbvrev.end[i]]*end^2
#rate[i] <- ynext[i]/yend[abbvrev.end[i]] - 1 
#}
}
"
