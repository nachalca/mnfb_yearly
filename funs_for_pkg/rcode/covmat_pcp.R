#' Parallel coordinate plot of a covariance matrix distribution
#' @usage covmat(x, xnames=NULL, colvar = 'dep')
#' @param x a list of covariance matrices
#' @param xnames a (optional) character vector with names for x
#' @param colorvar color variable of each covariance. Available options are effective dependence, 'dep', effective variance 'var', or an user privided index
#' @export
#' @examples
#' smpcovs <- rIW(n=100, d = 15, R = diag(4) )
#' covmat_pcp(x = smpcovs)
#'
covmat_pcp <- function(x, xnames=NULL, colvar = 'dep') {
  n <- length(x)
  k <- nrow(x[[1]])

  if (length(xnames) > 1) {
    names(x) <- xnames
  } else {
    names(x) <- 1:n
  }
  
  # compute variable to color lines
  effvar_fn <- function(x) svd( x )$d %>% log( ) %>% mean( ) %>% exp()
  effdep_fn <- function(x) 1 - svd( cov2cor(x) )$d %>% log( ) %>% mean( ) %>% exp()
  
  if (length(colvar) > 1) {
    cl.var <- scale(colvar) %>% as.numeric()
  }  else if (colvar == 'var') {
    cl.var <- sapply(x, effvar_fn )
  } else if (colvar == 'dep') {
    cl.var <- sapply(x, effdep_fn)
  }
  
  colvar.nm <- ifelse(length(colvar)>1, deparse(substitute(colvar) ) , colvar)
  color.dat <- data_frame(grp.var = names(x), cl.var = cl.var)
  
  # set up data for pc plot
  cov.data <- lapply(x, as.vector ) %>% bind_rows() %>% 
    mutate(row.n=rep(1:k, each=k), col.n=rep(1:k, times=k), 
           par.type = factor( row.n-col.n == 0, labels = c('cv', 'sg')   ), 
           param = paste( 'cov', rep(1:k, each=k), rep(1:k, times=k), sep='')) %>%
    select(-row.n, -col.n) %>%
    gather(grp.var, sg.val, -param, - par.type) %>%
    inner_join(color.dat)

  # produce final plot  
  cov.data %>% 
    ggplot( aes(param, sg.val, group=grp.var,fill=par.type) ) + 
    geom_point(size=2, shape=21, colour="grey50")  + 
    geom_line( aes( color = cl.var ) )  + 
    theme_bw() + 
    labs(x='', y='') + 
    theme(axis.text.x = element_text(angle = 60)) +
    scale_color_gradient2(name = colvar.nm, midpoint = median(cl.var) ) +
    scale_fill_manual( values=c("white", "black") ) + 
    guides( fill = FALSE ) 
}
