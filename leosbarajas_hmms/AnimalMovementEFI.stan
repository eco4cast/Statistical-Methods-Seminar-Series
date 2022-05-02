// This Stan program uses 'cmdstan' 

data {
  int<lower = 0> N;  // number of states
  int<lower = 1> T;  // number of observations
  int y[T]; // integer-valued
  
  //priors
  vector[N] lprior_shape;
  vector[N] lprior_scale;
}

parameters {
  simplex[N] theta[N];  // N x N tp
  ordered[N] lambda; // state-dependent parameters
  //vector<lower=0>[N] lambda; 
  simplex[N] init_dist;
}


transformed parameters {
  
  matrix[N, T] log_omega;
  matrix[N, N] Gamma;

  // build log_omega
  for (t in 1:T)
    for (n in 1:N) log_omega[n, t] = poisson_lpmf(y[t] | lambda[n]);

  // build Gamma
  for (n in 1:N) Gamma[n, ] = theta[n]';

}

model {
  
  // priors
  for(n in 1:N)
    lambda[n] ~ gamma(lprior_shape[n], lprior_scale[n]);
  
  target += hmm_marginal(log_omega, Gamma, init_dist);
  
}

generated quantities{

  matrix[N, T] states_prob = hmm_hidden_state_prob(log_omega, Gamma, init_dist);
  int states_pred[T] = hmm_latent_rng(log_omega, Gamma, init_dist);

   vector[T] y_pred;
   for (i in 1:T)
     y_pred[i] = poisson_rng(lambda[states_pred[i]]);

}

