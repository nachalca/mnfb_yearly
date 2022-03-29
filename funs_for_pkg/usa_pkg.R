
library(PriorCovmatrix)
? rSS

# small test
rIW(n=2, d = 6, R = diag(5) )

rSS(n=2, k = 5)
rSS(n=2, k = 5, prior_sg = 'ca')

rSS(n=2, k = 5, prior_cor = 'iw')
rSS(n=2, k = 5, prior_cor = 'iw', prior_sg = 'ca')

# 1 simulation to make some plots

sim1 <- rIW(n = 1e3, d = 6, R = diag(5) )


library(tidyverse)


# Data are created on mnfb repository, with the yearly_data.R code. 
bird.yeartotal <- read.csv('data/bird_yeartotal.csv') %>% 
  mutate( count.add = count*(count>0) + (count == 0)  ) %>% 
  mutate( ave.add  = ifelse(samples>0, (count.add/samples), 0)    , # compute average using the count.add (this step should be on yearly_data.R)
          yearc = year - 2000) 

# filter 10 species with more counts
lim <- bird.yeartotal %>% group_by( abbrev ) %>%
  summarise( total.count = sum(count) ) %>% 
  mutate(rr = rank(total.count)) %>%
  arrange(desc(rr)) %>% slice(11) %>%
  select(total.count) %>% as.numeric()

bird.top10 <- bird.yeartotal %>%
  group_by( abbrev ) %>%
  mutate( total.count = sum(count, na.rm = TRUE) ) %>%
  ungroup() %>% filter(total.count > lim) %>% select(-total.count)


bird.top10 %>% 
  mutate(abbrev = factor(abbrev)) %>%
  group_by(abbrev, year) %>%
  summarise( total.count = sum(count) ) %>%
  ggplot( aes (year, total.count, color = abbrev) ) + geom_point() + geom_line() 

smpcovs <- bird.top10 %>% 
  mutate(abbrev = factor(abbrev) ) %>%
  select(year, forest, abbrev, count) %>%
  spread(forest, count) %>% select(-abbrev) %>%
  nest(-year ) %>% 
  mutate( smp.cov = map( data, cov) )

smpcovs$smp.cov[[1]]
  
