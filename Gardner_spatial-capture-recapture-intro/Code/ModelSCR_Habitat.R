#############################################
##Statistical-Methods-Seminar-Series
##SCR webinar 
##Code to fit Model SCRhabitat to the bear data in NIMBLE
##sex effect on g0 and sigma
##Data from Gardner et al. 2009/2010
##dPointProcess by Nathan Hostetter
##Written by Beth Gardner 10/3/2026
#############################################

#Load packages
library(nimble)
library(coda)
library(MCMCvis)
library(terra)
library(raster)
library(reshape2)


########################### Nimble Function #################################
deregisterDistributions("dPointProcess")
#Make custom nimble function
dPointProcess <- nimbleFunction(run = function(x = double(1), pi=double(1), nPix=double(0), grid.dist = double(0),
                                                xlim=double(1), ylim=double(1), pixMat=double(2), log=integer(0, default = 0)) {
  returnType(double(0))
  if(x[1]>xlim[1] & x[1]<xlim[2] & x[2]>ylim[1] & x[2]<ylim[2])
  {#keep in state space
    #assign location to pixel
    sx <- trunc(((x[1]+abs(xlim[1]))/grid.dist)) + 1
    sy <- trunc(((x[2]+abs(ylim[1]))/grid.dist)) + 1
    pix <- pixMat[sx, sy]
    
    loglike <- dcat(pix, pi[1:nPix], log=TRUE)
  }else{
    loglike <- -Inf
  }
  if(log) return(loglike)
  else return(exp(loglike))
}
)

##############################################################################
##############################################################################
### SCR with forest cover as a density covariate in Nimble


#Define the nimble model
modelSCRhabitat<- nimbleCode({
  #priors
  g0~dbeta(1,1)  #detection probability
  sigma~dunif(0, 6)  #scale parameter, units of kms
  b0 ~ dnorm(-2, sd=2)
  b1 ~ dnorm(0, sd=2)
  
  mu[1:SS_size] <- exp(b0 + b1*habcov[1:SS_size])*pixelarea  #habitat model
  EN <- sum(mu[1:SS_size])   # expected number of individuals in state-space
  pi[1:SS_size] <- mu[1:SS_size]/EN
  psi <- EN / M #inclusion probability
  
  for(i in 1:M){       
         z[i] ~ dbern(psi) 
    #call new distribution to reduce computation time
    s[i, 1:2] ~ dPointProcess(pi=pi[1:SS_size], nPix=SS_size, grid.dist = grid.dist, 
                              xlim=xlim[1:2], ylim=ylim[1:2],
                              pixMat = pixMat[1:nPix_x, 1:nPix_y])
    
    d2[i,1:J]<- pow(s[i,1]-X[1:J,1],2) + pow(s[i,2]-X[1:J,2],2)
     g[i,1:J] <- z[i]*g0*exp(-d2[i, 1:J]/(2*sigma^2))  #detection function
    
    for(j in 1:J){
      y[i,j] ~ dbin(g[i,j], K) 
    }
  }
  
  N <- sum(z[1:M]) #abundance in the state space
})


##################################################
##Read in and prepare the data

#load in the bear data
load("Gardner_spatial-capture-recapture-intro/Data/beardata.rda")
yarray<-beardata$bearArray  #array of bear captures by trap and occasion
trapmat<-beardata$trapmat   #trap coordinates, in UTMs/1000 (units = kms)

nind<-dim(beardata$bearArray)[1]
K<-dim(beardata$bearArray)[3]
ntraps<-dim(beardata$bearArray)[2]

M=200 #M = observed + augmented individuals
nz<-M-nind  #nz is the number of augmented individuals

#create augmented array
Yaug <- array(0, dim=c(M,ntraps,K))
Yaug[1:nind,,]<-yarray
y<-apply(Yaug,1:2, sum)  #reduce to only M by ntraps, there are no time effects

#center the coordinates of the trap matrix
X=as.matrix(cbind((trapmat[,1]- mean(trapmat[,1])), (trapmat[,2]- mean(trapmat[,2]))))


# Load the forest cover data
forest_cover <- rast("Gardner_spatial-capture-recapture-intro/Data/forest_cover.tif")
pixelarea = prod(res(forest_cover))/(1000*1000)  ##dividing to make m to kms, should be .25km^2


grid.dist <- sqrt(pixelarea) # distance between centroids in kms
# number of rows and columns 
nPix_y  <- nrow(forest_cover)
nPix_x  <- ncol(forest_cover)
#making the xlim and ylim into kms
xlim <- (ext(forest_cover)[1:2] - mean(trapmat[,1])*1000)/1000
ylim <- (ext(forest_cover)[3:4] - mean(trapmat[,2])*1000)/1000
forest.df <- as.data.frame(forest_cover, xy = TRUE, na.rm = TRUE)


SS_size <- nrow(forest.df)
# assign pixel number as look-up style matrix
pixMat <- as.matrix(dcast(cbind(forest.df, 1:SS_size), x~y, fun.aggregate = sum, value.var ="1:SS_size"))[,-1]


#get mean activity centers for observed bears; create random initial values for augmented individuals
Sin<-matrix(NA, ncol=2, nrow=M)

for(i in 1:nind){
  Sin[i,]<- colMeans(y[i,]*X)
}
Sin[(nind+1):M,]<-cbind(runif(nz, xlim[1], xlim[2]), runif(nz, ylim[1], ylim[2]))



##################################################
##Set up data for NIMBLE

#Defining constants
constants <- list(M = M, 
                  K = K, #Occasions
                  J = ntraps, 
                  pixelarea = pixelarea,
                  SS_size = SS_size,
                  grid.dist=grid.dist, pixMat=pixMat,
                  nPix_x=nPix_x, nPix_y=nPix_y, xlim=xlim, ylim=ylim, X = X[,1:2])


#Define data
data <- list(y = y, habcov = forest.df[,3]/100)  #convert forest cover to decimal

#Set inits
inits <- list(sigma = runif(1,1,2),  #Inputs all these as initial values
                   g0=runif(1), 
                   s = Sin, 
                   z = c(rep(1,nind), rbinom(nz,1,0.5)),
                   b0 = -3, b1 = 0)

params<-c("N", "g0", "psi", "sigma", "b0", "b1")



#Fit the model
Rmodel <- nimbleModel(code = modelSCRhabitat,  #This builds the model in R, but is not compiled yet.
                           constants = constants, 
                           data = data, 
                           inits = inits)

conf <- configureMCMC(Rmodel,monitors=params, 
                           control = list(adaptInterval = 200), thin=1) 

#to adjust the samplers for b0 and b1 to improve convergence
#conf$removeSamplers(c("b0", "b1"))
#conf$addSampler(target = c("b0", "b1"), type = "RW_block")
#conf$addSampler(target = c("b0", "b1"), type = 'AF_slice')


Rmcmc <- buildMCMC(conf)#Building the chains
Cmodel <- compileNimble(Rmodel)#Compiles the model in c++
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)#Compile in c++

##This version takes about 1-5 minutes to run, will need to run longer to converge
start<-Sys.time()
samplesHabitat <- runMCMC(Cmcmc,           #Run compiled model 
                            niter = 5000,
                            nburnin = 1000,
                            nchains = 3)
end<-Sys.time()
end-start


## results, convergence checks
MCMCsummary(samplesHabitat, n.eff = TRUE, Rhat = TRUE)

MCMCtrace(object = samplesHabitat,
          pdf = FALSE,
          ind = TRUE)

#Make a plot of density
b0est<-MCMCsummary(samplesHabitat)[2,1]
b1est<-MCMCsummary(samplesHabitat)[3,1]

forest.df$z<-exp(b0est+b1est*forest.df[,3]/100)
forestests<-cbind(forest.df$x, forest.df$y, forest.df$z)
forestrast <- rasterFromXYZ(forestests)

plot(forestrast, col = terrain.colors(100))
