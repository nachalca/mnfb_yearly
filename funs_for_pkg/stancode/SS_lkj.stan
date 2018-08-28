// code to simulate from SS: cor = LKJ
data {
  int <lower=0> k; // covariance matrix dimension 
  real<lower=0> eta; // for LKJ prior
  real<lower=0> sigma_mu; 
  real<lower=0> sigma_sc;
  int prior_sigma ;
  }
  parameters {
    corr_matrix[k] Omega; // correlation matrix
    vector<lower=0>[k] sigma;
  }
model { 
  Omega ~ lkj_corr(eta );
  if (prior_sigma == 1)
    target += lognormal_lpdf(sigma | sigma_mu, sigma_sc);
  else if (prior_sigma == 2)
    target += cauchy_lpdf(sigma | sigma_mu, sigma_sc);
  }
generated quantities {
  cov_matrix[k] Sigma;
  Sigma = quad_form_diag(Omega, sigma);
}

