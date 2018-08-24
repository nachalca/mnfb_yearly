#' Random generation from IW distribution
#' @usage rIW(n, d, R)
#' @param n number of observations
#' @param d degrees of freedom IW parameter
#' @param R matrix IW parameter
#' @export
#' @example 
#' rIW(n=1e3, d = 3, R = diag(2) )

rIW <- function(n, d, R) { 
  rstan::sampling(rIW, data = list(d=d, R=R, k = ncol(R)), iter = n, chains = 1)
   }

