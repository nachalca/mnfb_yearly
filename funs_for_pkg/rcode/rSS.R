#' Random generation from IW distribution
#' @usage rIW(n, d, R)
#' @param n number of observations
#' @param k covariance matrix dimension 
#' @param prior_cor prior for correlation matrix
#' @param prior_sg  prior for marginal variances
#' @param eta degrees of freedom IW parameter
#' @param R matrix IW parameter
#' @param sigma_mu location of variances 
#' @param sigma_sc scale of variances
#' @export
#' @examples
#' rIW(n=4, d = 3, R = diag(2) )
#'

rSS <- function(n, k, prior_cor = 'lkj', prior_sg ='ln', 
                eta = k+1, R = diag(k), 
                sigma_mu=0, sigma_sc=1) {
  
  # obtain stan model object
  if (prior_cor == 'lkj') {
    mm <- get( stanmodels$SS_lkj )
  }  
  if (prior_cor == 'iw') {
    mm <- get( stanmodels$SS_iw)
  }
  
  # create data list
  dts <- list(k=k, sigma_mu=sigma_mu, eta=eta, sigma_sc=sigma_sc, 
              prior_sigma = 1*(prior_sg == 'ln') + 2*(prior_sg=='ca') )
  
  if (prior_cor == 'iw') dts$R = R
  
  # obtain samples
  out <- rstan::sampling(mm, data = dts,
                         iter = 2*n, chains = 1, refresh = 0) %>%
    rstan::extract(pars = 'Sigma')
  lapply(1:n, function(x) out$Sigma[x , , ] )
}

