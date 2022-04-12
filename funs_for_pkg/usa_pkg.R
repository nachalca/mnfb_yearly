
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

top10.cutpoint <- bird.yeartotal %>%
  group_by( abbrev ) %>%
  summarise( total.count = sum(count) ) %>% 
  arrange( desc(total.count) ) %>% 
  slice(11) %>% select(total.count) %>% as.numeric()
  
bird.top10 <- bird.yeartotal %>%
  group_by( abbrev ) %>%
  mutate( total.count = sum(count) ) %>% 
  filter(total.count > top10.cutpoint) %>% select(-total.count)

bird.top10 %>%
  group_by(abbrev, year) %>%
  summarise( total.count = sum(count) ) %>%
  ggplot(aes(year, total.count, color=abbrev)) + geom_point() + geom_line()

smpcovs <- bird.yeartotal %>%
  select(abbrev, forestN, year, ave.add) %>%
  spread(forestN, ave.add) %>% 
  group_by(abbrev) %>% nest(-year) %>%
  mutate( covmat = map( data, cov, use = 'complete.obs') )

covmat_pcp(x = smpcovs$covmat, xnames = smpcovs$abbrev)

zz <- bird.yeartotal %>% group_by(abbrev) %>%
  summarise( tt =  sum(count)) %>% pull(tt) %>%
  scale() %>% as.numeric()

covmat_pcp(x = smpcovs$covmat, xnames = smpcovs$abbrev, colvar =  zz)

covmat_pcp(x = smpcovs$covmat, xnames = smpcovs$abbrev, colvar =  'var')


smpcovs2 <- bird.yeartotal %>%
  select(abbrev, forestN, year, ave.add) %>%
  spread(forestN, ave.add) %>% filter(!is.na(Chequamegon) ) %>%
  group_by(year) %>% nest(-abbrev) %>%
  mutate( covmat = map( data, cov, use = 'complete.obs') )

covmat_pcp(x = smpcovs2$covmat, xnames = smpcovs2$year)

 
# plot of simulated inverse-wishart
pl <- expand.grid( d = c(5, 15), rho = c(0, .8) ) %>%
  apply( 1, function(x) {
    rr <- matrix(x[2], ncol=4, nrow =4 ) + diag(rep(1-x[2], 4))
    smpcovs <- rIW( n=100, d = x[1], R = rr )
    titu <- paste('d =', x[1], ',' , 'rho =', x[2], sep=' ')
    covmat_pcp(x = smpcovs) + ggtitle(titu) + guides(color = FALSE)
  })
grid.arrange(grobs = pl)




#-----
# a function to plot a coord par plot
library(GGally)
#-- args
x <- smpcovs$covmat
xnames <- as.character( smpcovs$abbrev )
cl.var <- bird.top10 %>% group_by(abbrev) %>%
  summarise( tt = as.numeric(sum(count)) )  %>%
  pull(tt) 
# ---

n <- length(x)
k <- nrow(x[[1]])
names(x) <- xnames

color.dat <- data_frame(grp.var = xnames, cl.var = scale(cl.var)%>% as.numeric())

cov.data <- lapply(x, as.vector ) %>% bind_rows() %>% 
  mutate(row.n=rep(1:k, each=k), col.n=rep(1:k, times=k), 
         par.type = factor( row.n-col.n == 0, labels = c('cv', 'sg')   ), 
         param = paste( 'cov', rep(1:k, each=k), rep(1:k, times=k), sep='')) %>%
    select(-row.n, -col.n) %>%
  gather(grp.var, sg.val, -param, - par.type) %>%
  inner_join(color.dat)

cov.data %>% 
  ggplot( aes(param, sg.val, group=grp.var,fill=par.type) ) + 
  geom_point(size=2, shape=21, colour="grey50")  + geom_line(aes(color=cl.var)) + 
  theme_bw() +
  scale_fill_manual(values=c("white", "black")) + 
  scale_color_gradient2()



