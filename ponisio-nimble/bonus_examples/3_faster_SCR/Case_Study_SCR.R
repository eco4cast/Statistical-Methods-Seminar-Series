## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
has_ggplot2 <- require(ggplot2)
has_mcmcplots <- require(mcmcplots)
has_coda <- require(coda)
has_nimbleSCR <- require(nimbleSCR)
generate_original_results <- TRUE
DataDir <- file.path("..","..","..","..","wolverine_data","DataScript") # Modify as needed. See below to get data


## ---- eval=FALSE--------------------------------------------------------------
## DataDir <- getwd() # Change this to control where the files will live
## download.file("https://datadryad.org/stash/downloads/file_stream/20780", file.path(DataDir, "wolverine_data.zip"))
## unzip(file.path(DataDir, "wolverine_data.zip"), "WolverineData.RData", exdir = DataDir)


## -----------------------------------------------------------------------------
load(file.path(DataDir, "WolverineData.RData"))

data <- list(y = my.jags.input$y,                           # data: individual x detector (grid cell) 
             z = my.jags.input$z,                           # indicators: 1 = real, NA = data augmented
             detector.xy = my.jags.input$detector.xy,       # grid of detector IDs
             habitat.mx = my.jags.input$habitat.mx,         # grid of habitat mask: 1 = valid habitat, 0 = not habitat
             ones = my.jags.input$OK,                       # dummy values for habitat mask
             lowerCoords = c(0,0),                          # coordinate ranges
             upperCoords = c(
               dim(my.jags.input$habitat.mx)[2],
               dim(my.jags.input$habitat.mx)[1]),
             trials = rep(1, dim(my.jags.input$detector.xy)[1]))
constants <- list(n.individuals = my.jags.input$n.individuals,
                  n.detectors = dim(my.jags.input$detector.xy)[1],
                  y.max = dim(my.jags.input$habitat.mx)[1],
                  x.max = dim(my.jags.input$habitat.mx)[2])
inits <- list(sxy = inits.1$sxy,                             # activity centers
              z = inits.1$z,
              p0 = 0.05,
              psi = 0.5,
              sigma = 6)


## ---- echo=FALSE--------------------------------------------------------------
# Do this here to get the figure generated.
set.seed(123)
DetectorIndex <- getLocalObjects(habitatMask = data$habitat.mx,
                                 coords = data$detector.xy,
                                 dmax = 38,
                                 resizeFactor = 24)

constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
constants$ResizeFactor <- DetectorIndex$resizeFactor
constants$n.cells <- dim(DetectorIndex$localIndices)[1]
constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
data$detectorIndex <- DetectorIndex$localIndices
data$nDetectors <- DetectorIndex$numLocalIndices
data$habitatIDDet <- DetectorIndex$habitatGrid


## -----------------------------------------------------------------------------
length(data$z)        # 1 for known individual, NA for data-augmented virtual individual
sum(data$z == 1, na.rm = TRUE)               # There are 196 known individuals.
dim(data$detector.xy) # There are 17266 "detectors", which are grid squares with (x, y) coordinates
dim(data$y)           # There are 874 x 17266 (individual x detector) data points
sum(data$y==1)        # The vast majority are 0s.  There are 372 detections.
dim(data$habitat.mx)  # Habitat mask: 1305 x 1089 grid points in a rectangle containing the study area.
sum(data$habitat.mx==1)  # 215796 grid squares of valid habitat in the study area with buffer


## ---- eval=FALSE--------------------------------------------------------------
## sink("SCR-LESS.jags") cat("model {
##     ##------------------------------------------------------------
##     ##----------    AC PLACEMENT      ---##
##     ##-----------------------------------##
##     for(i in 1:n.individuals){
##     sxy[i,1] ~ dunif(xy.bounds[i,1,1], xy.bounds[i,1,2])
##     sxy[i,2] ~ dunif(xy.bounds[i,2,1], xy.bounds[i,2,2])
##     pOK[i] <- habitat.mx[trunc(sxy[i,2])+1, trunc(sxy[i,1])+1]
##     OK[i] ~ dbern(pOK[i])
##     }#i
##     ##------------------------------------------------------------
##     ##----- DEMOGRAPHIC PROCESS ----##
##     ##------------------------------##
##     psi0 ~ dunif(0,1)
##     psi <- mean(psi1[])
##     for (i in 1:n.individuals){
##     psi1[i] <- 1-(1-psi0)^prop.habitat[i]
##     z[i] ~ dbern(psi1[i])
##     }#i
##     ##------------------------------------------------------------
##     ##---------- DETECTION PROCESS------##
##     ##----------------------------------##
##     p0 ~ dunif(0,1)
##     sigma ~ dunif(0,50)
##     alpha <- -1/(2*sigma*sigma)
##     #----- DETECTION PROCESS ------#
##     for (i in 1:n.individuals){
##     for (j in 1:n.detectors[i]){
##     d2[i,j] <- pow(sxy[i,1] - detector.xy[detector.index[i,j], 1] ,2) +
##                pow(sxy[i,2] - detector.xy[detector.index[i,j], 2], 2)
##     p[i,j] <- p0 * exp(alpha * d2[i,j])
##     y[i,detector.index[i,j]] ~ dbern(p[i,j]*z[i])
##     }#j
##     }#i
##     ##------------------------------------------------------------
##     ##---------- DERIVED PARAMETERS ----------##
##     ##----------------------------------------##
##     N <- sum(z[])
##     }",fill = TRUE)
## sink()


## -----------------------------------------------------------------------------
code <- nimbleCode({
  ## priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  p0 ~ dunif(0, 1)
  ## loop over individuals
  for(i in 1:n.individuals) {
    ## AC coordinates
    sxy[i,1] ~ dunif(0, x.max)
    sxy[i,2] ~ dunif(0, y.max)
    ## habitat constraint 
    ones[i] ~ dHabitatMask( s = sxy[i,1:2],
                            xmin = lowerCoords[1],
                            xmax = upperCoords[1],
                            ymin = lowerCoords[2],
                            ymax = upperCoords[2],
                            habitat = habitat.mx[1:y.max,1:x.max])
    ## latent virtual/real indicators
    z[i] ~ dbern(psi)
    ## likelihood
    y[i, 1:nMaxDetectors] ~ dbinomLocal_normal( detNums = nbDetections[i],
                                                detIndices = yDets[i,1:nMaxDetectors],
                                                size = trials[1:n.detectors],
                                                p0 = p0,
                                                s = sxy[i,1:2],
                                                sigma = sigma,
                                                trapCoords = detector.xy[1:n.detectors,1:2],
                                                localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                                localTrapsNum = nDetectors[1:n.cells],
                                                resizeFactor = ResizeFactor,
                                                habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                indicator = z[i])
  }
  ## derived quantity: total population size
  N <- sum(z[1:n.individuals])
})


## ---- eval=FALSE--------------------------------------------------------------
## DetectorIndex <- getLocalObjects(habitatMask = data$habitat.mx,
##                                  coords = data$detector.xy,
##                                  dmax = 38,
##                                  resizeFactor = 24) # This actually generates the map shown above.
## 
## constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]   # Coarse grid maximum values
## constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
## constants$ResizeFactor <- DetectorIndex$resizeFactor      # factor to scale ACs to coarse grid
## constants$n.cells <- dim(DetectorIndex$localIndices)[1]
## constants$maxNBDets <- DetectorIndex$numLocalIndicesMax   # Max number of local detectors relevant for a coarse grid cell
## data$detectorIndex <- DetectorIndex$localIndices          # Detector grid IDs
## data$nDetectors <- DetectorIndex$numLocalIndices          # Number of local detectors relevant for each coarse grid cell
## data$habitatIDDet <- DetectorIndex$habitatGrid            # Coarse grid IDs


## -----------------------------------------------------------------------------
ySparse <- getSparseY(x = my.jags.input$y)
data$y <- ySparse$y[,,1]                                # Number of detections only where detected
data$yDets <- ySparse$detIndices[,,1]                   # Indices where an animal was detected
data$nbDetections <- ySparse$detNums[,1]                # Number of detectors where an animal was detected
constants$nMaxDetectors <- ySparse$maxDetNums           # Maximum detections for any animal


## -----------------------------------------------------------------------------
dim(my.jags.input$y) # Size of riginal (dense) version of all detections
874 * 17266 / (2^20) # megabytes in original version
dim(data$y)          # Size of sparse version of all detections
874 * 8 / (2^20)     # megabytes in sparse version
8 / 17266            # Approx. reduction in iteration/memory work needed to use new version


## -----------------------------------------------------------------------------
## Consider animal 1
data$nbDetections[1]  # It was seen at 2 detectors
data$yDets[1,]        # At detectors 14941 and 14988
data$y[1,]            # Once at each detector (this data set has only 1s). "-1"s are fillers: do not use.
constants$nMaxDetectors  # Most detectors where one animal was seen was 8.


## -----------------------------------------------------------------------------
s <- c(324.6, 736.2)  # Example (x, y) coordinate for an activity center
# Activity center         --> coarse grid coordinates
# ResizeFactor scales from (x, y) to coarse grid
c(trunc(s[2]/constants$ResizeFactor) + 1, trunc(s[1]/constants$ResizeFactor) + 1)

# coarse grid coordinates --> coarse grid ID
data$habitatIDDet[trunc(s[2]/constants$ResizeFactor) + 1, trunc(s[1]/constants$ResizeFactor) + 1]
dim(data$habitatIDDet)                  # Matrix of coarse grid IDs for each (y,x)
range(data$habitatIDDet)                # Coarse grid IDs range from 0 to 509

# coarse grid ID          --> relevant detector indices (on detector grid, not coarse grid)
data$detectorIndex[168, 1:data$nDetectors[168] ] # The rest are zeros
dim(data$detectorIndex) 


## -----------------------------------------------------------------------------
with(constants, c(x.max, y.max)) # Be careful about (x, y) vs. (y, x)
s <- c(324.6, 736.2)             # Example (x, y) coordinate for an activity center
dim(data$habitat.mx)             # Habitat mask
data$habitat.mx[707, 416]        # How to look up whether an AC is in valid habitat, set up as (y, x)
data$habitat.mx[trunc(s[2]) + 1, trunc(s[1]) + 1] # How to do it from x = s[1], y = s[2]


## ---- eval=FALSE--------------------------------------------------------------
## dbinomLocal_normal # You can look at source code like any R function...


## ---- eval=FALSE--------------------------------------------------------------
## # ... but let's look at original source code with type declarations and comments
## dbinomLocal_normal <- nimbleFunction(
##   run = function( x = double(1),
##                   detNums = double(0, default = -999),
##                   detIndices = double(1),
##                   size = double(1),
##                   p0 = double(0, default = -999),
##                   p0Traps = double(1),
##                   sigma = double(0),
##                   s = double(1),
##                   trapCoords = double(2),
##                   localTrapsIndices = double(2),
##                   localTrapsNum = double(1),
##                   resizeFactor = double(0, default = 1),
##                   habitatGrid = double(2),
##                   indicator = double(0),
##                   lengthYCombined = double(0, default = 0),
##                   log = integer(0, default = 0)
##   ) {
##     ## Specify return type
##     returnType(double(0))
## 
##     ## Deal with cases where detection info is combined in one vector
##     if(detNums==-999){
##       detNums <- x[1]
##       nMaxDetectors <- (lengthYCombined-1)/2
##       detIndices1 <- x[(nMaxDetectors+2):lengthYCombined]
##       x1 <- x[2:(nMaxDetectors+1)]
##     }else{
##       x1 <- x
##       detIndices1 <- detIndices
##     }
## 
##     ## Shortcut if the current individual is not available for detection
##     if(indicator == 0){
##       if(detNums == 0){
##         if(log == 0) return(1.0)
##         else return(0.0)
##       } else {
##         if(log == 0) return(0.0)
##         else return(-Inf)
##       }
##     }
## 
##     ## Retrieve the index of the habitat cell where the current AC is
##     sID <- habitatGrid[trunc(s[2]/resizeFactor)+1, trunc(s[1]/resizeFactor)+1]
## 
##     ## Retrieve the indices of the local traps surrounding the selected habita grid cell
##     theseLocalTraps <- localTrapsIndices[sID,1:localTrapsNum[sID]]
## 
##     ## CHECK IF DETECTIONS ARE WITHIN THE LIST OF  LOCAL TRAPS
##     if(detNums > 0){
##       for(r in 1:detNums){
##         if(sum(detIndices1[r] == theseLocalTraps) == 0){
##           if(log == 0) return(0.0)
##           else return(-Inf)
##         }
##       }
##     }
## 
##     ## Calculate the log-probability of the vector of detections
##     alpha <- -1.0 / (2.0 * sigma * sigma)
##     logProb <- 0.0
##     detIndices1 <- c(detIndices1, 0)
##     count <- 1
## 
## 
##     if(p0==-999){# when p0 is provide through p0Traps
##       for(r in 1:localTrapsNum[sID]){
##         if(theseLocalTraps[r] == detIndices1[count]){
##           d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
##           p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
##           logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
##           count <- count + 1
##         }else{
##           d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
##           p <- p0Traps[theseLocalTraps[r]] * exp(alpha * d2)
##           logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
## 
##         }
##       }
##     }else{# when p0 is provide through p0
##       for(r in 1:localTrapsNum[sID]){
##         if(theseLocalTraps[r] == detIndices1[count]){
##           d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
##           p <- p0 * exp(alpha * d2)
##           logProb <-  logProb + dbinom(x1[count], prob = p, size = size[theseLocalTraps[r]], log = TRUE)
##           count <- count + 1
##         }else{
##           d2 <- pow(trapCoords[theseLocalTraps[r],1] - s[1], 2) + pow(trapCoords[theseLocalTraps[r],2] - s[2], 2)
##           p <- p0 * exp(alpha * d2)
##           logProb <- logProb + dbinom(0, prob = p, size = size[theseLocalTraps[r]], log = TRUE)
##         }
##       }
##     }
## 
## 
##     ## Return the probability of the vector of detections (or log-probability if required)
##     if(log)return(logProb)
##     return(exp(logProb))
##   })












## ---- echo=FALSE--------------------------------------------------------------
have_samples <- exists("samples")

