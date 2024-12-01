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
  for(i in 1:n) {
    target += 
    log(
      Phi(((y[i] - delta) * mu[i] - gamma) / sqrt((y[i] - delta) * sigma)) - 
      Phi(((y[i] - delta - 1) * mu[i] - gamma) / sqrt((y[i] - delta - 1) * sigma))
    );
  }
}
generated quantities {
  vector[m_upper - m_lower + 1] y_pred;
  for(i in m_lower:m_upper) {
    y_pred[i - m_lower + 1] = gamma / i + delta;
  }
}
