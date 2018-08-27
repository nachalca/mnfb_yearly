
# create package skeleton
# http://mc-stan.org/rstantools/articles/minimal-rstan-package.html
library(rstantools)
library(usethis)
library(roxygen2)


#1) create skeleton with rstan_package_skeleton

pkg.path <- '~/research/master_cc/pkgTest'
rstan_package_skeleton(path = pkg.path, 
                       stan_files = list.files('funs_for_pkg/stancode', full.names = T)
)

#2) Copy: 
#  - R files in R folder
#  - STAN files in src/stan_files folder

ll <- list.files('funs_for_pkg/rcode', full.names=T)

file.copy(ll, paste(pkg.path, '/R', sep = '') )


# 3) intall once using 
# (usando stan_files argument may not be necesary)
# https://github.com/r-lib/devtools/issues/1653
# https://github.com/stan-dev/rstanarm/issues/190
devtools::install(args = "--preclean")

#4) now can use
devtools::document()
devtools::install()
devtools::check()
