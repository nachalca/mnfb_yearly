data {
  int <lower=0> N;
  int <lower=0> k;  
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
}
parameters {
  vector[k] mu;
  cov_matrix[k] Sigma;
}
transformed parameters {
  real s1;
  real s2;
  real rho;
//  corr_matrix[k] Rho;
//  vector<lower=0>[k] sig;
  s1 = sqrt(Sigma[1,1]);
  s2 = sqrt(Sigma[2,2]);
  rho = Sigma[1,2]/(s1*s2);
// compute correlation matrix 
  // for (i in 1:k) sig[i] = 1/sqrt(Sigma[i,i]);
  // for (n in 1:k) {
  //   for (m in 1:n) {
  //     Rho[m,n] = sig[m] * sig[n] * Sigma[m,n]; 
  //   }
  // }
  // for (n in 1:k) {
  //    for (m in (n+1):k) {
  //       Rho[m,n] = Rho[n,m];
  //    }
  // } 

}
model {
  //for ( i in 1:k)  mu[i] ~ normal(0, 100);
  Sigma ~ inv_wishart(k+1, R);
  for (n in 1:N) y[n] ~ multi_normal(mu, Sigma);
}
