data {
  real<lower=0> eta;
  int<lower=1> nsims;
  int<lower=1> K;
}
parameters {
  corr_matrix[K] R;
}
model {
  R ~ lkj_corr(eta);
}
