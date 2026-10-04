#############################################
##Statistical-Methods-Seminar-Series
##SCR webinar 
##Code for running black bear data in secr
##Data from Gardner et al. 2009/2010
##Written by Beth Gardner updated 10/1/2026
#############################################

#load packages
library(secr)
library(terra)

#load in the bear data
load("Gardner_spatial-capture-recapture-intro/Data/beardata.rda")
traps<-read.csv("Gardner_spatial-capture-recapture-intro/secr/beartraps.csv")  
#beartraps is a trap file formated for secr and is in UTMs

##hair snares are a considered proximity detector in this example, 
##so we define that here 
trapdet<-read.traps(data=traps, detector = "proximity") 

##create a covariate to be sex
##beardata$flat uses the "session" to assign sex, here session 1 is female and 
##session 2 is male. We'll re-code this so session is 1 for all bears and create a
##covariate called "sex".

nb<-as.data.frame(beardata$flat[,1:4])
nb$sex<-beardata$flat[,1]
nb$Session <- rep(1, length(nb$Session))

bear.cap <- make.capthist(as.data.frame(nb), trapdet, covnames=c("sex"),
                          fmt = "trapID", noccasions = 8)

##Once the data are formatted, we can move forward to fitting models.
##Using a buffer of 20km here

## constant detection g0~1, sigma~1 and constant density
fit.0 <- secr.fit(bear.cap, hcov = "sex",buffer=20000, trace = FALSE,
                  model = list(D ~ 1, g0 ~ 1, sigma ~ 1))

#if you get an error that starts with:
#Error: Error in function ibeta_derivative<long double>
#(long double,long double,long double): Overflow Error10.
#add this "details" line into your secr.fit call:
"...details = list(fastproximity = FALSE),"

predict(fit.0)

## estimate of sex specific detection parameters and constant density
fit.h2 <- secr.fit(bear.cap, hcov = "sex", buffer=20000, trace = FALSE,
                     model = list(D ~ 1, g0 ~ h2, sigma ~ h2))
predict(fit.h2)


## estimate of learned behavior response in detection parameter and constant density
fit.b <- secr.fit(bear.cap, hcov = "sex", buffer=20000, trace = FALSE,
                   model = list(D ~ 1, g0 ~ b, sigma ~ 1))
predict(fit.b)


##Now we'll fit an inhomogeneous point process with percent forest cover

#Read in the forest cover data and create a habitat mask for secr
forest_cover<-rast("Gardner_spatial-capture-recapture-intro/Data/forest_cover.tif")

#convert spatial raster to a data frame for secr
forest.df <- as.data.frame(forest_cover, xy = TRUE, na.rm = TRUE)
forest.df$forest=forest.df$forest/100  #covert forest cover to a decimal

#convert the data frame into a secr mask object
habitat.mask <- read.mask(
  data = forest.df, 
  spacing = res(forest_cover)[1] # Match the cell size of your raster
)

##fit the model, specify density as a covariate for D

fit.habitat <- secr.fit(bear.cap, model = list(D ~ forest, g0 ~ 1, sigma ~ 1), 
                        mask = habitat.mask)

fit.habitat

#plot density surface
hold=predictDsurface(fit.habitat, mask = habitat.mask)
plot(hold)   
points(traps(bear.cap),col='black')

#Return the expected number of bears in the state space based on the provided
#habitat mask
region.N(fit.habitat,habitat.mask)

