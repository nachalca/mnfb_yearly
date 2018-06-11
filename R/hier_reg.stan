data {
  int<lower=0> N; // obs
  int<lower=1> K; // predictors
  int<lower=1> J; // groups
  int<lower=1,upper=J> jj[N];  // group for individual
  matrix[K,K] R; // for IW(J+1, R) prior
  vector[K] beta0; // for N(beta0, Sigma) prior
  matrix[N, K] x;
  vector[N] y;
}
parameters {
  cov_matrix[K] Sigma;
  vector[K] beta[J];
  real<lower=0> sigma;
}
model {
  Sigma ~ inv_wishart( K + 1, R );
  sigma ~ cauchy(0,1);

  for (j in 1:J) {
    beta[j] ~ multi_normal(beta0, Sigma);
  }
      
  for (n in 1:N) {
    y[n] ~ normal(x[n] * beta[jj[n]], sigma);
  }
}
