// code to simulate from IW(d, R) distribution, R is kxk matrix, d=k+1 implies unif on rhos
data {
  int <lower=0> d;
  matrix[k,k] R;
}
parameters {
  cov_matrix[k] Sig_iw;
}
model {
  Sig_iw ~ inv_wishart(d, R);   // IW prior 
}