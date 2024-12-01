data {
  int<lower=1> n;
  int<lower=1> m_lower;
  int<lower=m_lower> m_upper;
  matrix<lower=0>[n, 366] temp;
  vector<lower=0>[n] y;
}
parameters {
  real<lower=0> gamma;
  real<lower=0> sigma;
  real mu_a;
  real<lower=0> sigma_a;
  vector<lower=1, upper=80>[n] a;
}
model {
  vector[n] mu_sum;
  matrix[n, 366] w;
  
  for(i in 1:n) {
    for(j in 1:366) {
      if(j < y[i]) {
        w[i,j] = inv_logit(10 * (j - a[i]));
      } else {
        w[i,j] = 0;
      }
    }
    
  mu_sum[i] = dot_product(w[i], temp[i]);  
  target += 
    log(
      Phi((mu_sum[i] - gamma) / 
        sqrt((y[i] - a[i]) * sigma)) - 
      Phi((mu_sum[i] - temp[i, to_int(floor(y[i]))] - gamma) / 
        sqrt((y[i] - a[i] - 1) * sigma))
        );
  }
  a ~ normal(mu_a, sigma_a);
}
generated quantities {
  vector[m_upper - m_lower + 1] y_pred;
  for(i in m_lower:m_upper) {
    y_pred[i - m_lower + 1] = 90 * gamma / i + normal_rng(mu_a, sigma_a);
  }
}
