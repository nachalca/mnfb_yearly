#' Random generation from IW distribution
#' @usage rIW(n, d, R)
#' @param n number of observations
#' @param d degrees of freedom IW parameter
#' @param R matrix IW parameter
#' @export
#' @examples
#' rIW(n=4, d = 3, R = diag(2) )
#'
rIW <- function(n, d, R) {
  dts <- list(d=d, R=R, k = ncol(R))
  out <- rstan::sampling(stanmodels$rIW, data = dts,
                  iter = n, chains = 1, refresh = 0, warmup = 0,
                  algorithm = 'Fixed_param') %>%
    rstan::extract(pars = 'Sigma')
  lapply(1:n, function(x) out$Sigma[x , , ] )

   }

