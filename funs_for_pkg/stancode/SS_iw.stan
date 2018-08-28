// code to simulate from SS: cor = IW
data {
  int <lower=0> k; // covariance matrix dimension 
  real<lower=0> eta; // df for IW prior
  matrix[k, k] R; // location for IW prior
  real<lower=0> sigma_mu; 
  real<lower=0> sigma_sc;
  int prior_sigma ;
  }
  parameters {
    cov_matrix[k] xOmega; 
    vector<lower=0>[k] sigma;
  }
  transformed parameters {
  matrix[k , k] Omega; 
  vector[k] gamma2;
  vector[k] inv_gamma;
  gamma2    = diagonal(xOmega);
  inv_gamma =  inv_sqrt(gamma2) ;
  Omega     = quad_form_diag(xOmega, inv_gamma);
  }
 model { 
  xOmega  ~ inv_wishart(eta, R);
  if (prior_sigma == 1)
    target += lognormal_lpdf(sigma | sigma_mu, sigma_sc);
  else if (prior_sigma == 2)
    target += cauchy_lpdf(sigma | sigma_mu, sigma_sc);
  }
generated quantities {
  cov_matrix[k] Sigma;
  Sigma = quad_form_diag(Omega, sigma);
}

    

