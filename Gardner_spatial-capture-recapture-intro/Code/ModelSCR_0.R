#############################################
##Statistical-Methods-Seminar-Series
##SCR webinar 
##Model SCR0 for the bear data in NIMBLE
##Data from Gardner et al. 2009/2010
##Written by Beth Gardner updated 10/1/2026
#############################################


library(nimble)    #load nimble
library(coda)
library(MCMCvis)

#############################################
### SCR0 model in Nimble

modelSCR0<- nimbleCode( {

  #priors
psi~dunif(0,1) #inclusion probability
g0~dbeta(1,1)  #detection probability
sigma~dunif(0, 6)  #scale parameter, units of kms
sigma2<-sigma*sigma

for(i in 1:M){  #M observed + augmented individuals
 	z[i] ~ dbern(psi)     #indicator if individual is in the population
	s[i,1]~dunif(Xl,Xu)   #activity centers X coord
 	s[i,2]~dunif(Yl,Yu)   #activity centers Y coord

  for(j in 1:J){  #J traps
		d2[i,j]<- pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2)  #calculate the distance
		g[i,j]<- z[i]*g0*exp(-d2[i,j]/(2*sigma2))              #detection function
		y[i,j] ~ dbin(g[i,j], K) 
	}

}
N<-sum(z[])  #number of individuals in the state space
D<-N/area    #density of individuals in the state space
})

##################################################
##Read in and prepare the data

#load in the bear data
load("Data/beardata.rda")
yarray<-beardata$bearArray  #array of bear captures by trap and occasion
trapmat<-beardata$trapmat #trap coordinates, in UTMs/1000 (units = kms)

#yarray is a 3d array with 1:47 bears, 1:38 traps, and 1:8 occassions

nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=210 #M = observed + augmented individuals
nz<-M-nind  #nz is the number of augmented individuals

#create augmented array
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-yarray
y<-apply(Yaug,1:2, sum)  #reduce to only M by ntraps, there are no time effects

#center the coordinates of the trap matrix
X=as.matrix(cbind((trapmat[,1]- mean(trapmat[,1])), (trapmat[,2]- mean(trapmat[,2]))))

#set up the state-space; adding 8km to the edge of the trap array
Xl=min(X[,1]) - 8
Xu=max(X[,1]) + 8
Yl=min(X[,2]) - 8
Yu=max(X[,2]) + 8
areaX=(Xl-Xu)*(Yl-Yu)


#get mean activity centers for observed bears; create random initial values for augmented individuals
Sin<-matrix(NA, ncol=2, nrow=M)

for(i in 1:nind){
Sin[i,]<- colMeans(y[i,]*X)
}
Sin[(nind+1):M,]<-cbind(runif(nz, Xl, Xu), runif(nz,Yl,Yu))

##################################################
##Set up data for NIMBLE

data<-list(y=y)
constants <- list(M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
params<-c('psi','g0','sigma', 'N', 'D')

inits = list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, 
		sigma=runif(1,2,3),g0=runif(1))

##Using the nimbleMCMC function here to consolidate nimble code.
##This version takes about 1-5 minutes to run
start<-Sys.time()
samples <- nimbleMCMC(
    code = modelSCR0,  
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
MCMCsummary(samples, n.eff = TRUE, Rhat = TRUE)

MCMCtrace(object = samples,
          pdf = FALSE,
          ind = TRUE)


