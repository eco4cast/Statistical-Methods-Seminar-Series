# 19 April 2023
# iSDM that combines presence-only data with replicated counts
# data simulation & code taken from Ch. 10.6 of Applied Hierarchical Modeling in Ecology, Vol 2
# using data simulation function built into the AHMbook R package

library(AHMbook)
library(nimble)
library(MCMCvis)
library(raster)

set.seed(111)

# simDataDK - from AHM book package
dat <- AHMbook::simDataDK(
  sqrt.npix = 100, 
  alpha = c(-2, -1), 
  beta = c(6, 0.5), 
  drop.out.prop.pb = 0, 
  quadrat.size = 5,
  gamma = c(0, -1.5), 
  nquadrats = 25, 
  nsurveys = 3, 
  show.plot = TRUE )

logarea <- log(dat$s.area / dat$npix)
y <- numeric(dat$npix)
y[dat$pixel.id.det] <- 1

data <- list(
  y = y,
  logarea = logarea, 
  xcov = dat$xcov, 
  wcov = dat$wcov,
  C = dat$countData[,5:7], 
  covarX = dat$countData[,"x"],
  covarW = dat$countData[,"w"],
  area = dat$s.area * (dat$quadrat.size^2 / dat$npix),
  nquadrats = ncell( dat$squad ))

constants <- list(
  npix = dat$npix,
  nsite = dat$nquadrats, 
  nsurveys = dat$nsurveys)

code <- nimble::nimbleCode({
  
  beta0 ~ dnorm(0, sd = 2)
  beta1 ~ dnorm(0, sd = 2)
  
  alpha0 ~ dnorm(0, sd = 2)
  alpha1 ~ dnorm(0, sd = 2)
  
  # Submodel 1 for presence-only data
  #-----------------------------------
  
  for( i in 1:npix ) {
    y[i] ~ dbern( psi [i] )
    cloglog( psi[i] ) <- logarea + log(lambda[i]) + log(p1[i])
    log( lambda[i] ) <- beta0 + beta1 * xcov[i]
    logit( p1[i] ) <- alpha0 + alpha1 * wcov[i]
  }
  
  # Submodel 2 for replicated counts
  #---------------------------------
  gamma0 <- logit( mean.p )
  mean.p ~ dunif(0, 1)
  gamma1 ~ dnorm(0, sd = 2)
  
  for( i in 1:nsite) {
    N[i] ~ dpois( area * lam[i] )
    log( lam[i] ) <- beta0 + beta1 * covarX[i]
    for( j in 1:nsurveys ) {
      C[i, j] ~ dbin( p2[i], N[i] )
    }
    logit( p2[i] ) <- gamma0 + gamma1 * covarW[i]
  }
  
  # Calculate abundance as a derived variable
  Nppm <- sum( lambda[1:npix] * exp(logarea) )
  Nnmix <- (sum( N[1:nsite]) / nsite ) * nquadrats
})

Nst <- apply(data$C, 1, max) + 1
inits <- function(){
  list(
    N = Nst, 
    beta0 = rnorm(1, 0, 1), 
    beta1 = rnorm(1, 0, 1), 
    alpha0 = rnorm(1, 0, 1), 
    alpha1 = rnorm(1, 0, 1), 
    mean.p = runif(1, 0.1, 0.9), 
    gamma1 = rnorm(1, 0, 1)
  )
}

params <- c("alpha0", "alpha1", "beta0", "beta1", "gamma0", "gamma1", "mean.p", "Nppm", "Nnmix")

# MCMC settings
nc <- 1
nburn <- 2000  
ni <- nburn + 3000
nt <- 3

# takes ~2 minutes to run on my machine
start <- Sys.time()
# run the model!
out <- nimble::nimbleMCMC(
  code = code, 
  constants = constants, 
  data = data, 
  inits = inits(),
  monitors = params, 
  niter = ni, 
  nburnin = nburn, 
  nchains = nc, 
  thin = nt )
end <- Sys.time()
end - start

MCMCvis::MCMCsummary(out)