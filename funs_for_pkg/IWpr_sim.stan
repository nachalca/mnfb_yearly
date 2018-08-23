// code to simulate from IW(d, R) distribution, R is kxk matrix, d=k+1 implies unif on rhos
data {
  int <lower=0> nu;
  int <lower=0> d;
  matrix[d,d] R;
}
parameters {
  cov_matrix[d] Sig_iw;
}
model {
  Sig_iw ~ inv_wishart(nu, R);   // IW prior 
}
