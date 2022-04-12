data {
  int<lower=0> N; // number of obs
  int<lower=1> K; // number of predictor variables
  int<lower=1> J; // groups
  int<lower=1,upper=J> jdx[N];  // group for individual
  matrix[N, K] x;
  vector[N] y;
}
parameters {
  corr_matrix[K] Omega;
  vector<lower=0>[K] sigma_b;
  vector[K] beta[J];
  real<lower=0> sigma;
}
model {
  Omega ~ lkj_corr( K );
  sigma_b ~ cauchy(0,1);
  sigma ~ cauchy(0,1);

  for (j in 1:J) {
    beta[j] ~ multi_normal(0, quad_form_diag(Omega, sigma_b));
  }
      
  for (n in 1:N) {
    y[n] ~ normal(x[n] * beta[jdx[n]], sigma);
  }
}
