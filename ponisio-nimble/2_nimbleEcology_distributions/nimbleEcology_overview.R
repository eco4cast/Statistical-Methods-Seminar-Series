## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
library(nimble)
library(coda)
has_nimbleEcology <- require(nimbleEcology)
has_compareMCMCs <- require(compareMCMCs)
if(!has_nimbleEcology)
  message("This module will use nimbleEcology, which you don't have installed.")
doDHMMexample <- TRUE


## ---- echo=FALSE--------------------------------------------------------------
if(!has_compareMCMCs) message("To run the code as given, you'll need to install compareMCMCs (see workshop installation information).  Otherwise you can run code via normal nimble workflows.")


## ---- eval=FALSE--------------------------------------------------------------
## y[i, 1:T] ~ dOcc_v(probOcc = psi[i],
##                    probDetect = p[i, 1:T], len = T)


## ---- eval=FALSE--------------------------------------------------------------
## ## Not indexing transition and observation probabilities by individual, for simplicity of example
## for(t in 1:T) {
##   for(i in 1:num_states) {
##     for(j in 1:num_states) {
##       T[i, j, t] <- some_calculation() # Involving survival and transition probabilities
##                                        # that may depend on explanatory variables or random effects
##     }
##   }
## }
## ## Some similar definition of O, observation probabilities
## for(k in 1:num_individuals) {
##   for(t in 2:T) {
##     z[k, t] ~ dcat(T[z[t-1], 1:num_states, t])
##   }
##   for(t in 1:T) {
##     y[k, t] ~ dcat(O[z[t], 1:num_states, t])
##   }
## }


## -----------------------------------------------------------------------------
orchids_example_dir <- file.path("..", "..", "content", "examples","orchids") # modify for your own system
CH <- as.matrix(read.table(file.path(orchids_example_dir, "orchids.txt"), sep=" ", header = FALSE))
head(CH) ## Glimpse of data
n_occasions <- dim(CH)[2]
## Compute vector with occasion of first capture
f <- numeric()
for (i in 1:dim(CH)[1])
  f[i] <- min(which(CH[i,]!=0))

## Modification for NIMBLE:
CH <- CH[which(f!=11), ]  ## remove all individuals not seen until the last occasion: They contribute no information and create problems for the dDHMM version.
## Reset f from the reduced CH
f <- numeric()
for (i in 1:dim(CH)[1])
  f[i] <- min(which(CH[i,]!=0))

## Recode CH matrix: note, a 0 is not allowed by WinBUGS!
## 1 = seen vegetative, 2 = seen flowering, 3 = not seen 
rCH <- CH  # Recoded CH 
rCH[rCH==0] <- 3 

## Function to create known latent states z 
known_state_ms <- function(ms, notseen){
    ## notseen: label for not seen 
    state <- ms 
    state[state==notseen] <- NA 
    for (i in 1:dim(ms)[1]){
        m <- min(which(!is.na(state[i,]))) 
        state[i,m] <- NA 
    }
    return(state) 
}
## Bundle data 
orchids_data <- list(y = rCH,
                  f = f,
                  n_occasions = dim(rCH)[2],
                  nind = dim(rCH)[1],
                  z = known_state_ms(rCH, 3)) 

ms_init_z <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   states <- max(ch, na.rm = TRUE)
   known.states <- 1:(states-1)
   v <- which(ch==states)
   ch[-v] <- NA
   ch[v] <- sample(known.states, length(v), replace = TRUE)
   return(ch)
   }

## Initial values 
orchids_inits <- function(){
    list(s = runif((dim(rCH)[2]-1), 0, 1),
         z = ms_init_z(rCH, f)
         )}



## -----------------------------------------------------------------------------
orchids_code <- nimbleCode({
  ## -------------------------------------------------
  ## Parameters:
  ## s: survival probability 
  ## psiV: transitions from vegetative 
  ## psiF: transitions from flowering 
  ## psiD: transitions from dormant 
  ## -------------------------------------------------
  ## States (S):
  ## 1 vegetative 
  ## 2 flowering 
  ## 3 dormant 
  ## 4 dead 
  ## Observations (O):  
  ## 1 seen vegetative 
  ## 2 seen flowering 
  ## 3 not seen 
  ## -------------------------------------------------
  ## Priors and constraints 
  ## Survival: uniform 
  for (t in 1:(n_occasions-1)){  
    s[t] ~ dunif(0, 1) 
  }
  ## Transitions: gamma priors 
  for (i in 1:3){
    a[i] ~ dgamma(1, 1) 
    psiD[i] <- a[i]/sum(a[1:3]) 
    b[i] ~ dgamma(1, 1) 
    psiV[i] <- b[i]/sum(b[1:3]) 
    c[i] ~ dgamma(1, 1) 
    psiF[i] <- c[i]/sum(c[1:3]) 
  }
  ## Define state-transition and observation matrices 	
  for (i in 1:nind){
    ## Define probabilities of state S(t+1) given S(t) 
    for (t in 1:(n_occasions-1)){
      ps[1,i,t,1] <- s[t] * psiV[1]
      ps[1,i,t,2] <- s[t] * psiV[2]
      ps[1,i,t,3] <- s[t] * psiV[3]
      ps[1,i,t,4] <- 1-s[t]
      ps[2,i,t,1] <- s[t] * psiF[1]
      ps[2,i,t,2] <- s[t] * psiF[2]
      ps[2,i,t,3] <- s[t] * psiF[3]
      ps[2,i,t,4] <- 1-s[t]
      ps[3,i,t,1] <- s[t] * psiD[1]
      ps[3,i,t,2] <- s[t] * psiD[2]
      ps[3,i,t,3] <- s[t] * psiD[3]
      ps[3,i,t,4] <- 1-s[t]
      ps[4,i,t,1] <- 0 
      ps[4,i,t,2] <- 0 
      ps[4,i,t,3] <- 0 
      ps[4,i,t,4] <- 1 
      ## Define probabilities of O(t) given S(t) 
      po[1,i,t,1] <- 1 
      po[1,i,t,2] <- 0 
      po[1,i,t,3] <- 0 
      po[2,i,t,1] <- 0 
      po[2,i,t,2] <- 1 
      po[2,i,t,3] <- 0 
      po[3,i,t,1] <- 0 
      po[3,i,t,2] <- 0 
      po[3,i,t,3] <- 1 
      po[4,i,t,1] <- 0 
      po[4,i,t,2] <- 0 
      po[4,i,t,3] <- 1 
    } #t 
  } #i 
  ## Likelihood
  for (i in 1:nind){
    ## Define latent state at first capture 
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n_occasions){
      ## State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4]) 
      ## Observation process: draw O(t) given S(t) 
      y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:3]) 
    } #t 
  } #i
})
set.seed(123)
orchids_info <- list(code=orchids_code, constants=orchids_data, inits=orchids_inits())




## -----------------------------------------------------------------------------
orchids_code_faster <- nimbleCode({
  ## -------------------------------------------------
  ## Parameters:
  ## s: survival probability 
  ## psiV: transitions from vegetative 
  ## psiF: transitions from flowering 
  ## psiD: transitions from dormant 
  ## -------------------------------------------------
  ## States (S):
  ## 1 vegetative 
  ## 2 flowering 
  ## 3 dormant 
  ## 4 dead 
  ## Observations (O):  
  ## 1 seen vegetative 
  ## 2 seen flowering 
  ## 3 not seen 
  ## -------------------------------------------------
  ## Priors and constraints 
  ## Survival: uniform 
  for (t in 1:(n_occasions-1)){  
    s[t] ~ dunif(0, 1) 
  }
  ## Transitions: gamma priors 
  for (i in 1:3){
    a[i] ~ dgamma(1, 1) 
    psiD[i] <- a[i]/sum(a[1:3]) 
    b[i] ~ dgamma(1, 1) 
    psiV[i] <- b[i]/sum(b[1:3]) 
    c[i] ~ dgamma(1, 1) 
    psiF[i] <- c[i]/sum(c[1:3]) 
  }
  ## Define state-transition and observation matrices 	
  ## Define probabilities of state S(t+1) given S(t) 
  for (t in 1:(n_occasions-1)){
    ps[1,t,1] <- s[t] * psiV[1]
    ps[1,t,2] <- s[t] * psiV[2]
    ps[1,t,3] <- s[t] * psiV[3]
    ps[1,t,4] <- 1-s[t]
    ps[2,t,1] <- s[t] * psiF[1]
    ps[2,t,2] <- s[t] * psiF[2]
    ps[2,t,3] <- s[t] * psiF[3]
    ps[2,t,4] <- 1-s[t]
    ps[3,t,1] <- s[t] * psiD[1]
    ps[3,t,2] <- s[t] * psiD[2]
    ps[3,t,3] <- s[t] * psiD[3]
    ps[3,t,4] <- 1-s[t]
    ps[4,t,1] <- 0 
    ps[4,t,2] <- 0 
    ps[4,t,3] <- 0 
    ps[4,t,4] <- 1 
  }
  ## Define probabilities of O(t) given S(t) 
  po[1,1] <- 1 
  po[1,2] <- 0 
  po[1,3] <- 0 
  po[2,1] <- 0 
  po[2,2] <- 1 
  po[2,3] <- 0 
  po[3,1] <- 0 
  po[3,2] <- 0 
  po[3,3] <- 1 
  po[4,1] <- 0 
  po[4,2] <- 0 
  po[4,3] <- 1 
  ## Likelihood
  for (i in 1:nind){
    ## Define latent state at first capture 
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):n_occasions){
      ## State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], t-1, 1:4]) 
      ## Observation process: draw O(t) given S(t) 
      y[i,t] ~ dcat(po[z[i,t], 1:3]) 
    } #t 
  } #i
})






## -----------------------------------------------------------------------------
orchids_code_DHMM <- quote({
  ## -------------------------------------------------
  ## Parameters:
  ## s: survival probability 
  ## psiV: transitions from vegetative 
  ## psiF: transitions from flowering 
  ## psiD: transitions from dormant 
  ## -------------------------------------------------
  ## States (S):
  ## 1 vegetative 
  ## 2 flowering 
  ## 3 dormant 
  ## 4 dead 
  ## Observations (O):  
  ## 1 seen vegetative 
  ## 2 seen flowering 
  ## 3 not seen 
  ## -------------------------------------------------
  ## Priors and constraints 
  ## Survival: uniform 
  for (t in 1:(k-1)){  
    s[t] ~ dunif(0, 1)
  }
  ## Transitions: gamma priors 
  for (i in 1:3){
    a[i] ~ dgamma(1, 1) 
    psiD[i] <- a[i]/sum(a[1:3]) 
    b[i] ~ dgamma(1, 1) 
    psiV[i] <- b[i]/sum(b[1:3]) 
    c[i] ~ dgamma(1, 1) 
    psiF[i] <- c[i]/sum(c[1:3]) 
  }
  ## Define state-transition and observation matrices 	
  for (t in 1:(k-1)) {
    T[1,1,t] <- s[t] * psiV[1]
    T[1,2,t] <- s[t] * psiV[2]
    T[1,3,t] <- s[t] * psiV[3]
    T[1,4,t] <- 1-s[t]
    T[2,1,t] <- s[t] * psiF[1]
    T[2,2,t] <- s[t] * psiF[2]
    T[2,3,t] <- s[t] * psiF[3]
    T[2,4,t] <- 1-s[t]
    T[3,1,t] <- s[t] * psiD[1]
    T[3,2,t] <- s[t] * psiD[2]
    T[3,3,t] <- s[t] * psiD[3]
    T[3,4,t] <- 1-s[t]
    T[4,1,t] <- 0
    T[4,2,t] <- 0
    T[4,3,t] <- 0
    T[4,4,t] <- 1
  }
  O[1,1] <- 1 
  O[1,2] <- 0 
  O[1,3] <- 0 
  O[2,1] <- 0 
  O[2,2] <- 1 
  O[2,3] <- 0 
  O[3,1] <- 0 
  O[3,2] <- 0 
  O[3,3] <- 1 
  O[4,1] <- 0 
  O[4,2] <- 0 
  O[4,3] <- 1 
  for(i in 1:nind) {
    for(j in 1:4)
      init[i, j] <- y2[i, f[i]]==j ## y2 is the same as y, to avoid a cycle in graph (i.e. y[i, f[i]] depends on init[i, 1:3] which depends on y[i, f[i]])
  }
  for (i in 1:nind) {
    y[i, f[i]:k] ~ dDHMM(init = init[i, 1:4],
                         probObs = O[1:4,1:3],
                         probTrans = T[1:4,1:4,f[i]:(k-1)],  ## Note: possible length 1 issue
                         len = k-f[i]+1,
                         checkRowSums = 0)
  }
})
orchids_constants_DHMM <- list(f=f,
                               k=dim(rCH)[2],
                               nind=dim(rCH)[1])
orchids_data_DHMM <- list(y = rCH, y2 = rCH)
orchids_inits_DHMM <- function() {
  list(s = runif((dim(rCH)[2]-1), 0, 1),
       a = rep(1,3),
       b = rep(1,3),
       c = rep(1,3))
}
set.seed(123)
orchids_info_DHMM <- list(code=orchids_code_DHMM, 
                          constants=orchids_constants_DHMM,
                          data = orchids_data_DHMM,
                          inits=orchids_inits_DHMM())






## ---- eval=FALSE--------------------------------------------------------------
## # Code for your own results:
## make_MCMC_comparison_pages(c(orchids_result, orchids_result_faster, orchids_result_DHMM),
##                            dir = "your_orchid_results",
##                            modelName = "orchids")

