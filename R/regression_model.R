
# Develop a motivation example for isec 2018
# hierarchical regression model

library(tidyverse); library(rstan)

# Data are created on mnfb repository, with the yearly_data.R code. 
bird.yeartotal <- read.csv('data/bird_yeartotal.csv') %>%
  mutate( ave.add  = count.add/samples, # compute average using the count.add (this step should be on yearly_data.R)
          yearc = year - 2000) 

################
# species info  FALTA !!!!!!
code <- read.csv('~/Documents/mnfb/data/nrri_bird_code.csv',header=T)
#######################

# Overall trend within forest 
bird.yeartotal <- bird.yeartotal %>% 
  group_by(abbrev) %>% mutate( sp.total = sum(count)) %>%
  ungroup() %>% mutate(sp.abun = sp.total > quantile(sp.total, prob = .7) ) 

bird.yeartotal %>%
  group_by(year, forestN) %>% mutate( bird.total = sum(count)) %>%
  ggplot(aes(x=year, y=bird.total, shape=forestN, color=forestN) ) + 
  geom_point() + geom_line()

bird.yeartotal %>%
  ggplot() + geom_point(aes(yearc, yearc^2))

obs.cfs <- bird.yeartotal %>% 
  group_by(abbrev) %>% nest() %>% 
  mutate( mds = map(data, function(d) lm( log(count.add) ~ yearc + I(yearc^2), data=d)  ) ) %>%
  mutate(cfs = map(mds, 
                   function(m) { data.frame( coef(m) ) %>% 
                     rownames_to_column(var='nm')  %>%
                       spread(nm, coef.m.)
                     }
                    ) )

obs.cfs$cfs[[1]]
obs.cfs %>% unnest(cfs) %>% with(cor(`I(yearc^2)`, yearc))

obs.cfs %>% unnest(cfs) %>% 
  ggplot() + geom_point(aes(`I(yearc^2)`, yearc), size=I(3) ) + 
  labs(x='quadratic term', y='linear term') + 
  theme(axis.title.x =  element_text(size=I(20) ), axis.title.y = element_text(size=I(20)) ) +
  ggsave( filename = 'report/figs/ols_regs.pdf', height = 7, width = 7) 

# Bayesian hierarchical model
# with rstan?? 
reg.iw <- stan_model(file = 'R/hier_reg.stan')
reg.lkj <- stan_model(file = 'R/hier_reg_lkj.stan')

dt.ls <- with(bird.yeartotal, list(
  N=nrow(bird.yeartotal), K = 3, J=length(unique(abbrev)), 
  jj = as.integer( factor(abbrev) ),
  R = diag(3), beta0 = c(0,0,0), 
  x = cbind(1, yearc, yearc^2) %>% as.matrix(),
  y = log(count.add)
))

res.reg.iw <- sampling(reg.iw, data = dt.ls)  
saveRDS(res.reg.iw, 'isec2018/reg_example.rds')

plot(res.reg.iw, plotfun='rhat')
plot(res.reg.iw, plotfun='ess')

rhos <- extract(res.reg.iw, pars='Sigma')$Sigma %>%
  apply( 1, function(z) z[2,3]/(sqrt(z[2,2]*z[3,3]))  ) %>%
  as_data_frame() %>% set_names(nm = 'smp')
  

ggplot(rhos) + geom_density(aes(x = smp))

ggplot(rhos) + geom_histogram(aes(x=smp, y=..density..)) +
labs(x='correlation', y = '') + 
  theme_bw() + 
  theme(axis.title.x =  element_text(size=I(20)), axis.text.x =  element_text(size=I(20))) +
  ggsave( filename = 'report/figs/reg_rho_hist.pdf', height = 7, width = 7) 


res.reg.lkj <- sampling(reg.lkj, data = dt.ls)  
saveRDS(res.reg.lkj, 'isec2018/reg_example_lkj.rds')


plot(res.reg.lkj, plotfun = 'rhat')

rhos.lkj <- extract(res.reg.lkj, pars='Omega')$Omega %>%
  apply( 1, function(z) z[2,3]  ) %>%
  as_data_frame() %>% set_names(nm = 'smp')


ggplot(rhos.lkj) + geom_density(aes(x = smp))



