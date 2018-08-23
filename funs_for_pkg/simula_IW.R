#' Funtion to do something
#' @usage simula_IW(n, d, R)
#' @param n dfjasdlkf
#' @param d kjdsklfas
#' @param R something
#' @export
#' @example 
#' simula_IW(n=1e3, d=3, R=diag(2) )

simula_IW <- function(n, d, R) { 
  rstan::sampling(mod, data = list(d=d, R=R, k = ncol(R)), iter = n, chains = 1)
 # mod <- stan_model(file = 'funs_for_pkg/IWpr_sim.stan')
   }

