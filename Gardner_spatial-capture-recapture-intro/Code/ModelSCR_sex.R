#############################################
##Statistical-Methods-Seminar-Series
##SCR webinar 
##Code to fit Model SCRsex to the bear data in NIMBLE
##sex effect on g0 and sigma
##Data from Gardner et al. 2009/2010
##Written by Beth Gardner updated 10/1/2026
#############################################


library(nimble)   
library(coda)
library(MCMCvis)

##############################################################################
### SCR with sex as det. covariate model in Nimble

modelSCRsex<- nimbleCode( {
psi~dunif(0,1) #inclusion probability
pi~dunif(0,1)  #sex ratio (probability of being female)

for(t in 1:2){
g0[t]~dbeta(1,1)  #detection probability
sigma[t]~dunif(0, 6)  
sigma2[t]<-sigma[t]*sigma[t]
}

for(i in 1:M){
   z[i] ~ dbern(psi)   #indicator if individual is in the population
 SEX[i] ~ dbern(pi)    #0/1 indicator if individual is male or female
SEX2[i] <- SEX[i] + 1  #add 1 so that we can index by the individual sex below
 s[i,1] ~ dunif(Xl,Xu) #activity centers X coord
 s[i,2] ~ dunif(Yl,Yu) #activity centers Y coord

for(j in 1:J){
d2[i,j]<- pow(s[i,1]-X[j,1],2) + pow(s[i,2]-X[j,2],2)
 g[i,j]<- z[i]*g0[SEX2[i]]*exp(-d2[i,j]/(2*sigma2[SEX2[i]]))
 y[i,j] ~ dbin(g[i,j],K)
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
trapmat<-beardata$trapmat   #trap coordinates, in UTMs/1000 (units = kms)
sex<-beardata$sex           #sex of individual bears

nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3 ]
ntraps<-dim(beardata$bearArray)[2]

M=220
nz<-M-nind

#create augmented array
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-yarray
y<-apply(Yaug,1:2, sum)  #reduce to only M by ntraps

#center the coordinates of the trap matrix
X=as.matrix(cbind((trapmat[,1]- mean(trapmat[,1])), (trapmat[,2]- mean(trapmat[,2]))))

#set up the state-space
Xl=min(X[,1]) - 8
Xu=max(X[,1]) + 8
Yl=min(X[,2]) - 8
Yu=max(X[,2]) + 8
areaX=(Xl-Xu)*(Yl-Yu)

#get mean activity centers for observed bears; create initial values for remaining s
Sin<-matrix(NA, ncol=2, nrow=M)

for(i in 1:nind){
Sin[i,]<- colMeans(y[i,]*X)
}
Sin[(nind+1):M,]<-cbind(runif(nz, Xl, Xu), runif(nz,Yl,Yu))

#create vector with the indicator of sex, all NAs for unobserved bears
SEX<-c(sex-1, rep(NA, nz))

#create an initialization vector for sex that has NA for observed bears and 0/1 for unobserved bears
SEXin=c(rep(NA, nind), rbinom(nz, 1,0.5))

##################################################
##Set up data for NIMBLE

data<-list(y=y,SEX=SEX)
constants <- list(M=M,K=K, J=ntraps, Xl=Xl, Yl=Yl, Xu=Xu, Yu=Yu, X=X, area=areaX)
params<-c('psi','g0','N', 'D', 'sigma', 'pi')
inits = list(z=c(rep(1,nind), rbinom(nz,1,0.5)),psi=runif(1), s=Sin, SEX=SEXin,
		pi=runif(1), sigma=runif(2,2,3),g0=runif(2))

##This version takes about >5 minutes to run, definitely needs longer, but it's not terrible!
start<-Sys.time()

samplesSex <- nimbleMCMC(
    code = modelSCRsex,  
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
MCMCsummary(samplesSex, n.eff = TRUE, Rhat = TRUE)

MCMCtrace(object = samplesSex,
          pdf = FALSE,
          ind = TRUE)




