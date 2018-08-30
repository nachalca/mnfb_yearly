
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
bird.top10 <- bird.yeartotal %>%
  group_by( nrricode ) %>%
  mutate( total.count = sum(count), rr = rank(total.count, ties.method = 'random') ) %>%
  filter(rr < 11)  %>% select(-total.count, -rr)


bird.top10 %>%
  group_by(nrricode, year) %>%
  summarise( total.count = sum(count) ) %>%
  with( table(year, nrricode) )



