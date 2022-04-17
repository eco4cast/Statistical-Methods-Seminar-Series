library(nimble)
ZIPcode <- nimbleCode({
  p ~ dunif(0,1)
  lambda ~ dunif(0,10)
  for (i in 1:N)
    y[i] ~ dZIP(lambda, zeroProb = p) ## Note NIMBLE allows R-like named-parameter syntax
})
dZIP <- nimbleFunction(
 run = function(x = double(), lambda = double(), # Type declarations
                zeroProb = double(), log = integer(0, default = 0)) {
   returnType(double())                           # return type declaration
   prob_x <- zeroProb*(x==0) + (1-zeroProb) * dpois(x, lambda, log = FALSE)
   if (log) return(log(prob_x))
   return(prob_x)
 })
set.seed(1)
y <- rbinom(50, size = 1, prob = 0.7) * rpois(50, lambda = 3)
mZIP <- nimbleModel(ZIPcode, data = list(y = y), constants = list(N = 50))
mcmcZIP <- buildMCMC(mZIP)
CmZIP <- compileNimble(mZIP) ## SOME PROBLEM WITH rZIP here.
samples <- nimbleMCMC(ZIPcode, data = list(y = y), 
                      constants = list(N = 50), niter = 1000)

