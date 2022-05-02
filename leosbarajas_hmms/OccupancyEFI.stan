data {
    int<lower=0> Tlen; // length of the time series
    int y[Tlen]; // observations
    int<lower=1> N; // number of states
    int index[Tlen];
}

parameters {
    real<lower = 0, upper = 1> p; // detection probability

    // regression coefficients for transition probabilities
    simplex[N] tpm[N];
    simplex[N] init_dist; 
}  

model {
    vector[N] logp;
    vector[N] logptemp;
    vector[N] log_gamma_tr[N];
    
    // priors
    p ~ beta(10, 1);

    // transpose
    for(i in 1:N)
        for(j in 1:N)
            log_gamma_tr[j,i] = log(tpm[i,j]);

    // likelihood computation
    logp[1] = log(init_dist[1]) + log(1-y[1]);
    logp[2] = log(init_dist[2]) + bernoulli_lpmf(y[1]| p);
    
    for (t in 1:Tlen) {
      
      if( t==1 || index[t] != index[t-1]){
        logp[1] = log(init_dist[1]) + log(1-y[t]);
        logp[2] = log(init_dist[2]) + bernoulli_lpmf(y[t]| p);
      } else {
        for(n in 1:N){
          logptemp[n] = log_sum_exp(to_vector(log_gamma_tr[n]) + logp);
        }
        logptemp[1] = logptemp[1] + log(1-y[t]);
        logptemp[2] = logptemp[2] + bernoulli_lpmf(y[t]|p);
     
        logp = logptemp;
      }
      
      if(t == Tlen || index[t+1]!= index[t])
          target += log_sum_exp(logp);
    }
}

