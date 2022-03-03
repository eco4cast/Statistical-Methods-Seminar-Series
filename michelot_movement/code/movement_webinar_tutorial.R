
# Plotting package
library(ggplot2)
theme_set(theme_bw())
# Movement modelling packages
library(momentuHMM)
library(foieGras)
library(adehabitatLT)
# GIS packages
library(sf)
library(sp)
# Load functions for later
source("utility_functions.R")

# Load data from Movebank URL
URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/10255/",
              "move.981/ThermochronTracking%20Elephants%20Kruger%202007.csv")
raw <- read.csv(url(URL))

# Keep relevant columns: ID, time, longitude, latitude, temperature
data_all <- raw[, c(9, 3, 4, 5, 6)]
colnames(data_all) <- c("ID", "time", "lon", "lat", "temp")
data_all$time <- as.POSIXct(data_all$time, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

# Just keep 2000 observations to save time with model fitting
data <- data_all[c(which(data_all$ID == unique(data_all$ID)[5])[1:1000],
                   which(data_all$ID == unique(data_all$ID)[6])[1:1000]),]

head(data)

# Plot latitude vs longitude (Mercator projection)
ggplot(data, aes(lon, lat, col = ID)) +
    geom_point(size = 0.5) + geom_path() +
    coord_map("mercator")

# Longitude vs time
ggplot(data, aes(time, lon, col = ID)) +
    geom_point(size = 0.5) + geom_path()

# Latitude vs time
ggplot(data, aes(time, lat, col = ID)) +
    geom_point(size = 0.5) + geom_path()

# Project to UTM
llcoord <- st_as_sf(data[, c("lon", "lat")], coords = c("lon", "lat"), 
                    crs = CRS("+proj=longlat +datum=WGS84"))
utmcoord <- st_transform(llcoord, crs = CRS("+proj=utm +zone=35 +datum=WGS84"))

# Add Easting-Northing to data (in km)
data[, c("x", "y")] <- st_coordinates(utmcoord)/1000

# Plot Northing vs Easting
ggplot(data, aes(x, y, col = ID)) + 
    geom_point(size = 0.5) + geom_path() +
    coord_equal()

# Table of time intervals in data
plot(table(diff(data$time)), xlim = c(0, 300), 
     xlab = "time interval (min)", ylab = "count")

# Use function from utility_function.R to split data at gaps > 2 hours
data_split <- split_at_gap(data = data, max_gap = 2*60, shortest_track = 24*60)
ggplot(data_split, aes(x, y, col = ID)) + 
    geom_point(size = 0.5) + geom_path() +
    coord_equal()

# Create adehabitat trajectory padded with NAs
data_ade <- setNA(ltraj = as.ltraj(xy = data_split[, c("x", "y")], 
                                   date = data_split$time, 
                                   id = data_split$ID), 
                  date.ref = data_split$time[1], 
                  dt = 30, tol = 5, units = "min")

# Transform back to dataframe
data_na <- ld(data_ade)[, c("id", "x", "y", "date")]
colnames(data_na) <- c("ID", "x", "y", "time")

# Add temperatures for non-missing locations
data_na$temp <- NA
data_na$temp[which(!is.na(data_na$x))] <- data_split$temp

head(data_na, 10)

# Prepare data for HMM (compute step lengths and turning angles)
data_hmm1 <-prepData(data_na, type = "UTM", covNames = "temp")

head(data_hmm1, 10)

# Observation distributions (step lengths and turning angles)
dist <- list(step = "gamma", angle = "vm")

# Initial parameters
# (step mean 1, step mean 2, step SD 1, step SD 2) and (angle concentration 1, angle concentration 2)
Par0_2s <- list(step = c(0.05, 0.2, 0.05, 0.2), angle = c(0.1, 3))

# Fit a 2-state HMM
hmm1 <- fitHMM(data_hmm1, nbStates = 2, dist = dist, Par0 = Par0_2s)

# Print parameter estimates
hmm1

# Plot estimated distributions and state-coloured tracks
plot(hmm1, breaks = 25, ask = FALSE)

# Initial parameters for 3-state model
Par0_3s <- list(step = c(0.02, 0.1, 0.3, 0.02, 0.1, 0.3), 
                angle = c(0.01, 0.1, 3))

# Fit 3-state HMM
hmm2 <- fitHMM(data_hmm1, nbStates = 3, dist = dist, Par0 = Par0_3s)

hmm2

plot(hmm2, breaks = 25, ask = FALSE)

# Get most likely sequence of states (Viterbi algorithm)
head(viterbi(hmm2))

# Save most likely state sequences from 2-state and 3-state models
data_hmm1$state_2st <- factor(viterbi(hmm1))
data_hmm1$state_3st <- factor(viterbi(hmm2))

# Plot tracks, coloured by states
ggplot(data_hmm1, aes(x, y, col = state_2st, group = ID)) +
    geom_point(size = 0.5) + geom_path() +
    coord_equal()
ggplot(data_hmm1, aes(x, y, col = state_3st, group = ID)) +
    geom_point(size = 0.5) + geom_path() +
    coord_equal()

# Fit 2-state HMM with temperature covariate (linear or quadratic effect)
hmm3 <- fitHMM(data_hmm1, nbStates = 2, dist = dist, 
               Par0 = Par0_2s, formula = ~temp)
hmm4 <- fitHMM(data_hmm1, nbStates = 2, dist = dist, 
               Par0 = Par0_2s, formula = ~temp+I(temp^2))

# Compare models using AIC
AIC(hmm3, hmm4)

# Plot estimated distributions and transition probabilities as functions of temperature
plot(hmm3, plotTracks = FALSE, ask = FALSE, plotCI = TRUE)

# Plot stationary state probabilities as functions of temperature
plotStationary(hmm3, plotCI = TRUE)

# Plot pseudo-residuals for 2-state and 3-state models
plotPR(hmm1)
plotPR(hmm2)

# Change data to format expected by foieGras
data_foieGras <- data_split[,c("ID", "time", "lon", "lat")]
colnames(data_foieGras)[1:2] <- c("id", "date")
# Add column for location quality class (G for "GPS")
data_foieGras$lc <- "G"
# Change order of columns as expected by foieGras
data_foieGras <- data_foieGras[,c(1, 2, 5, 3, 4)] 

# Fit state-space model to predict regular locations on 0.5h grid
ssm <- fit_ssm(d = data_foieGras, time.step = 0.5)

# Data frame of regularised locations
pred <- grab(ssm, what = "predicted", as_sf = FALSE)
data_reg <- as.data.frame(pred[, 1:4])
colnames(data_reg)[1:2] <- c("ID", "time")

plot(ssm, ask = FALSE)

# Get step lengths and turning angles from regularised data
data_hmm2 <- prepData(data = data_reg, type = "LL", coordNames = c("lon", "lat"))

head(data_hmm2, 10)

# Fit 2-state HMM to regularised data
hmm4 <- fitHMM(data_hmm2, nbStates = 2, dist = dist, Par0 = Par0_2s)

plot(hmm4, ask = FALSE)

# Predict locations on 30-min grid using crawl (through momentuHMM wrapper)
crw_out <- crawlWrap(obsData = data_split, timeStep = "30 min", 
                     Time.name = "time", coord = c("x", "y"))
data_hmm3 <- prepData(data = crw_out)

# Fit 2-state HMM to regularised data
hmm5 <- fitHMM(data_hmm3, nbStates = 2, dist = dist, Par0 = Par0_2s)

# Fit HMM using multiple imputation
hmm6 <- MIfitHMM(miData = crw_out, nSims = 10, ncores = 3, nbStates = 2, 
                 dist = dist, Par0 = Par0_2s)
