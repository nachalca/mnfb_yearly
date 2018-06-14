
# Develop a motivation example for isec 2018
# hierarchical regression model

library(tidyverse); library(rstan)
#library(rlist)

# Data are created on mnfb repository, with the yearly_data.R code. 
bird.yeartotal <- read.csv('data/bird_yeartotal.csv') %>% 
  mutate( count.add = count*(count>0) + (count == 0)  ) %>% 
  mutate( ave.add  = ifelse(samples>0, (count.add/samples), 0)    , # compute average using the count.add (this step should be on yearly_data.R)
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


# Function to compute separate models of 3 responses
sep_md <- function(df, resp) {
  df %>%
    mutate( response = resp) %>%
    group_by(abbrev) %>% nest() %>% 
    mutate( mds = map(data, function(d) lm( response ~ yearc + I(yearc^2), data=d)  ) ) %>%
    mutate(cfs = map(mds, 
                     function(m) { data.frame( coef(m) ) %>% 
                         rownames_to_column(var='nm')  %>%
                         spread(nm, coef.m.)
                     }
    ) )
}

# Function to compute an hierarchical model with IW prior
hier_md <- function( df, resp, mm ) {
  dt.ls <-df %>% mutate( response = resp) %>% 
    with(list(
      N=nrow(df), K = 3, J=length(unique(abbrev)), 
      jj = as.integer( factor(abbrev) ),
      R = diag(3), beta0 = c(0,0,0), 
      x = cbind(1, yearc, yearc^2) %>% as.matrix(), 
      y = response
    ))
  sampling(mm, data = dt.ls)  
}


res_sep_sum <- sep_md(df = bird.yeartotal, resp = bird.yeartotal$count.add)
res_sep_log <- sep_md(df = bird.yeartotal, resp = log(bird.yeartotal$count.add) )
res_sep_ave <- sep_md(df = bird.yeartotal, resp = bird.yeartotal$ave.add)

res_sep <- bind_rows(sum = res_sep_sum, log=res_sep_log, ave=res_sep_ave, .id = 'resp')

rm(res_sep_ave, res_sep_log, res_sep_sum)

res_sep %>% unnest(cfs) %>%
  ggplot() + geom_point(aes(`I(yearc^2)`, yearc), size=I(3) ) + 
    facet_wrap(~resp, scales = 'free') + 
  labs(x='quadratic term', y='linear term') + 
  theme(axis.title.x =  element_text(size=I(20) ), axis.title.y = element_text(size=I(20)) )
  
  
  ggsave( filename = 'report/figs/ols_regs.pdf', height = 7, width = 7) 

res_sep %>% unnest(cfs) %>% group_by(resp) %>%
  summarise( rr = cor( `I(yearc^2)`, yearc) )







##########################################################################
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

dt.ls <- with(bird.yeartotal, list(
  N=nrow(bird.yeartotal), K = 3, J=length(unique(abbrev)), 
  jj = as.integer( factor(abbrev) ),
  R = diag(3), beta0 = c(0,0,0), 
  x = cbind(1, yearc, yearc^2) %>% as.matrix()
))

# Use 2 models: using IW and sep strategy based on lkj 
reg.iw <- stan_model(file = 'R/hier_reg.stan')
reg.lkj <- stan_model(file = 'R/hier_reg_lkj.stan')


res.reg.iwsum <- sampling(reg.iw, data = dt.ls %>% list.append(y=bird.yeartotal$count) )  
saveRDS(res.reg.iwsum, 'isec2018/reg_iwsum.rds')

res.reg.iwlog <- sampling(reg.iw, data = dt.ls %>% list.append(y=log(bird.yeartotal$count)) )  
saveRDS(res.reg.iwlog, 'isec2018/reg_iwlog.rds')


res.reg.lkj <- sampling(reg.lkj, data = dt.ls)  
saveRDS(res.reg.lkj, 'isec2018/reg_example_lkj.rds')

plot(res.reg.iw, plotfun='rhat')
plot(res.reg.iw, plotfun='ess')

rhos <- extract(res.reg.iw, pars='Sigma')$Sigma %>%
  apply( 1, function(z) z[2,3]/(sqrt(z[2,2]*z[3,3]))  ) %>%
  as_data_frame() %>% set_names(nm = 'smp')

ggplot(rhos) + geom_histogram(aes(x=smp, y=..density..)) +
labs(x='correlation', y = '') + 
  theme_bw() + 
  theme(axis.title.x =  element_text(size=I(20)), axis.text.x =  element_text(size=I(20))) +
  ggsave( filename = 'report/figs/reg_rho_hist.pdf', height = 7, width = 7) 

rr <- extract(res.reg.lkj, pars='Omega')$Omega %>%
  apply( 1, function(z) z[2,3]  ) %>%
  as_data_frame() %>% set_names(nm = 'smp') %>%
  bind_rows( rhos, .id = 'met' ) %>%
  mutate(met = factor(met, labels = c('lkj', 'iw') ))

table(rr$met)
ggplot(rr) + geom_density(aes(x = smp, color = met))  +
  labs(x='correlation', y = '') + 
  theme_bw() + 
  theme(axis.title.x =  element_text(size=I(20)), axis.text.x =  element_text(size=I(20))) +
  ggsave( filename = 'report/figs/reg_rho_den.pdf', height = 7, width = 7) 

