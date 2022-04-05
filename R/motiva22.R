# Develop a motivation example for paper 2022
# 1 correlation coef between 2 species in 1 forest

library(tidyverse); library(rstan)
rstan_options(auto_write = TRUE)

# Data are created on mnfb repository, with the yearly_data.R code. 
bird.yeartotal <- read.csv('data/bird_yeartotal.csv') %>% 
  mutate( count.add = count*(count>0) + (count == 0)  ) %>% 
  mutate( ave.add  = ifelse(samples>0, (count.add/samples), 0)    , # compute average using the count.add (this step should be on yearly_data.R)
          yearc = year - 2000) 

################
# species info  FALTA !!!!!!
# code <- read.csv('~/Documents/mnfb/data/nrri_bird_code.csv',header=T)
#######################

with(bird.yeartotal, table(abbrev))
with(bird.yeartotal, tapply(count.add, abbrev, mean)) |> sort()


# usando 2 species
subset(bird.yeartotal, subset = forest == 9030) |> with(tapply(count.add, abbrev, mean)) |> sort()

dd <- subset(bird.yeartotal, 
             subset = forest == 9030,
             select=c('abbrev','ave.add','year', 'forest') ) |> 
  pivot_wider(names_from = abbrev, values_from = ave.add)

cor(dd, dd$WTSP) |> data.frame() |> setNames(nm='rr') |> subset(subset=rr > .7)


#We choose: VEER, WTSP in 9030 forest
dd2sp <- subset(bird.yeartotal, 
                subset = abbrev %in% c('VEER', 'WTSP') & forest == 9030,
                select=c('abbrev','ave.add', 'count.add', 'year', 'forest') ) |> 
  pivot_wider(names_from = abbrev, values_from = c(count.add, ave.add) ) 

ggplot(dd2sp) + geom_point(aes(count.add_VEER, count.add_WTSP, color=factor(forest)))

K <- 2
res.tot <- stan(file = 'R/normal_iw.stan', 
            data = list(y = dd2sp[,c('count.add_WTSP', 'count.add_VEER')], 
                        N = nrow(dd2sp), R = diag(K), k=K, mu0 = rep(0,K))
            )

# sqrt(var(dd2sp[,c('count.add_WTSP', 'count.add_VEER')]))
# cor(dd2sp[,c('count.add_WTSP', 'count.add_VEER')])
# res.tot
# plot(res.tot, pars='rho')


res.ave <- stan(file = 'R/normal_iw.stan', 
                data = list(y = dd2sp[,c('ave.add_WTSP', 'ave.add_VEER')], N = nrow(dd2sp), R = diag(K), k=K, mu0 = rep(0,K)), 
                pars = c('s1', 's2', 'rho'))

# var(dd2sp[,c('ave.add_WTSP', 'ave.add_VEER')])
# cor(dd2sp[,c('ave.add_WTSP', 'ave.add_VEER')])
# res.ave
# plot(res.ave, pars='rho')

true.rho <- data.frame( 
  cor(dd2sp[,c('count.add_WTSP', 'count.add_VEER')])[1,2], 
  cor(dd2sp[,c('ave.add_WTSP', 'ave.add_VEER')] )[1,2] ) |> 
  setNames(nm = c('Total', 'Mean')) |>  
  pivot_longer( names_to = 'response', values_to = 'rho', cols=1:2 )

p2 <- data.frame( 
  ave = extract(res.ave, pars='rho'), 
  cnt = extract(res.tot, pars='rho')
) |> 
  setNames(nm = c('Mean', 'Total')) |>  
  pivot_longer( names_to = 'response', values_to = 'rho', cols=1:2 ) |> 
  ggplot() + 
    geom_histogram(aes(x=rho), bins = 50) + facet_grid(response~., scale='free_y') + 
    geom_vline(data=true.rho, aes(xintercept=rho), linetype='dashed', size=I(1), color='red') +
  labs(x = 'Correlation coefficient', y='') +
  theme_bw()
  
p1 <- subset(bird.yeartotal, 
                subset = abbrev %in% c('VEER', 'WTSP') & forest == 9030,
                select=c('abbrev','ave.add', 'count.add', 'year') ) |> 
  setNames(nm = c('abbrev', 'mean', 'total', 'year')) |> 
  pivot_longer(cols = 2:3, names_to = 'response', values_to = 'yy') |> 
  pivot_wider(names_from = abbrev, values_from = yy ) |> 
  ggplot() + 
  geom_point(aes(WTSP, VEER)) + facet_wrap(response~., scale='free', ncol = 1)

library(patchwork)
p1 + p2

ggsave( filename = 'figs/motiva22.pdf', height = 7, width = 7) 

#=================================================
# load('data/simdata.Rdata')
# m_iw  <- stan_model(model_code=sim.iw)
# simdata <- subset(data, ns==size)                       
# mod_iw <-  dlply(simdata[simdata$ms =='iw', ], .(sim,r,s,ns),
#                  runstan.sim, prm=prms, .parallel=parallel)                        
# ms=c('iw', 'siw', 'ht', 'ss')
# #Run simulations for Bivariate case
# data2 <- data.frame( ms=rep(ms, each=nrow(simdata.2)),
#                      rbind(simdata.2,
#                            simdata.2,
#                            simdata.2,simdata.2) )



