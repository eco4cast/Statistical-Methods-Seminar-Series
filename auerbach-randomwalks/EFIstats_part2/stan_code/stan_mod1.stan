data {
  int<lower=1> n;
  int<lower=1> m_lower;
  int<lower=m_lower> m_upper;
  vector<lower=0>[n] mu;
  vector<lower=0>[n] y;
}
parameters {
  real<lower=0> gamma;
  real<lower=0> sigma;
  real delta;
}
model {
  y ~ normal(gamma / mu + delta, gamma * pow(sigma, 2) / pow(mu, 3));
}
generated quantities {
  vector[m_upper - m_lower + 1] y_pred;
  for(i in m_lower:m_upper) {
    y_pred[i - m_lower + 1] = gamma / i + delta;
  }
}
