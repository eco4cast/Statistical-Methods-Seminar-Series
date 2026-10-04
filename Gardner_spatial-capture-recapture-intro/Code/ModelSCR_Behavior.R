#############################################
##Statistical-Methods-Seminar-Series
##SCR webinar 
##Code to fit Model SCRB to the bear data in NIMBLE
##behavior effect by time, but not trap
##Data from Gardner et al. 2009/2010
##Written by Beth Gardner updated 10/1/2026
#############################################


library(nimble)    #load nimble
library(coda)
library(MCMCvis)

##############################################################################
### SCR with behavioral covariate model in Nimble
### After initial capture anywhere in the array, different baseline detection 
### probability is estimate. 

modelSCRB<- nimbleCode( {
  
alpha0~dnorm(0,.1)
alpha1~dnorm(0,.1)

sigma~dunif(0, 4)
sigma2<-sigma*sigma
psi~dbeta(1,1)

for(i in 1:M){
 z[i] ~ dbern(psi)
 s[i,1]~dunif(Xl,Xu)
 s[i,2]~dunif(Yl,Yu)

 d2[i,1:J]<- pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2)
 
for(k in 1:K){
  logit(g0[i,k])<- alpha0 + alpha1*C[i,k]
  g[i,1:J,k]<- z[i]*g0[i,k]*exp(-d2[i,1:J]/(2*sigma2))
  
for(j in 1:J){
  y[i,j,k] ~ dbin(g[i,j,k],1)
  }
 }

}
N<-sum(z[])
D<-N/area
})


##################################################
##Read in and prepare the data

#load in the bear data
load("Gardner_spatial-capture-recapture-intro/Data/beardata.rda")
yarray<-beardata$bearArray  #array of bear captures by trap and occasion
trapmat<-beardata$trapmat #trap coordinates, in UTMs/1000

nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3 ]
ntraps<-dim(beardata$bearArray)[2]

M=300
nz<-M-nind

#create augmented array
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-yarray
y<-apply(Yaug,1:2, sum)

#center the coordinates of the trap matrix
X=as.matrix(cbind((trapmat[,1]- mean(trapmat[,1])), (trapmat[,2]- mean(trapmat[,2]))))

#set up the state-space

Xl=min(X[,1]) - 8
Xu=max(X[,1]) + 8
Yl=min(X[,2]) - 8
Yu=max(X[,2]) + 8
areaX=(Xl-Xu)*(Yl-Yu)

###Create matrix to indicate if a bear has been captured previously at any trap

C=matrix(0, M, K)
for(k in 1:7){
    C[rowSums(Yaug[,,k]) >0, (k+1):8] <- 1
}


#get mean activity centers for observed bears; create initial values for remaining s
Sin<-matrix(NA, ncol=2, nrow=M)

for(i in 1:nind){
  Sin[i,]<- colMeans(y[i,]*X)
}
Sin[(nind+1):M,]<-cbind(runif(nz, Xl, Xu), runif(nz,Yl,Yu))


##################################################
##Set up data for NIMBLE

data<-list(y=Yaug)
constants <- list(M=M,K=K,C=C, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
params<-c('psi','alpha0','alpha1','N', 'D', 'sigma')

inits = list(z=c(rep(1,nind), rbinom(nz,1,0.5)), psi=runif(1), s=Sin,
             sigma=runif(1,2,3),alpha0=runif(1), alpha1=runif(1))


##This version takes about 15 minutes to run, does not converge, will need longer run
start<-Sys.time()
samplesB <- nimbleMCMC(
  code = modelSCRB,  
  data=data,
  constants = constants, 
  inits = inits,
  monitors = params,
  niter = 4000,    
  nburnin = 1000,
  nchains = 3,
  samplesAsCodaMCMC = TRUE,
  thin = 1)
end<-Sys.time()


end-start

## results, convergence checks
MCMCsummary(samplesB, n.eff = TRUE, Rhat = TRUE)

MCMCtrace(object = samplesB,
          pdf = FALSE,
          ind = TRUE)


