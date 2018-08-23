// code to simulate from IW(d, R) distribution, R is kxk matrix, d=k+1 implies unif on rhos
data {
  int <lower=0> k;
  int <lower=k> d;
  matrix[k,k] R;
}
model { }
generated quantities {
  cov_matrix[k] Sig_iw;
  Sig_iw = inv_wishart_rng(d, R);   // IW prior 
}
