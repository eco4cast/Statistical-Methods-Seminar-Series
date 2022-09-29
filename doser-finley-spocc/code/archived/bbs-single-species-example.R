# TODO: there is some bug going on here that you need to figure out. 
# bbs-single-species-example.R: this script fits a single-species occupancy model 
#                               using data from the North American Breeding Bird
#                               Survey. 
# Data source citation:   
#    Pardieck, K., Ziolkowski, D., Jr., Lutmerding, M., Aponte, V., & Hudson, M.- A. (2020). 
#    North american breeding bird survey dataset 1966– 2019. U.S. Geological Survey 
#    data release. https://doi.org/10.5066/ P9J6QUF6
rm(list = ls())
library(spOccupancy)
# For plotting and summarizing results
library(MCMCvis)
library(tidyverse)
library(sf)
library(stars)
library(pals)
# Set working directory if necessary
# setwd()
# Set seed for same results
set.seed(23)

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.veery)
load("data/bbs2019Veery.rda")
str(data.veery)

# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, single-species occupancy model with elevation 
out.1 <- PGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2), 
               det.formula = ~ scale(day) + (1 | obs), 
               data = data.veery, 
               n.samples = 5000, 
               n.thin = 4, 
               n.burn = 3000, 
               n.chains = 3,
               n.report = 500)
summary(out.1)
# Fit a non-spatial, single-species occupancy model with elevation and
# an intercept to allow occupancy to vary across Bird Conservation Regions.
out.2 <- PGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + factor(bcr) - 1, 
               det.formula = ~ scale(day) + (1 | obs), 
               data = data.veery, 
               n.samples = 5000, 
               n.thin = 4, 
               n.burn = 3000, 
               n.chains = 3,
               n.report = 500)
summary(out.2)

# Fit a spatial, single-species occupancy model with elevation and 
# a spatially-varying intercept using a NNGP and 10 nearest neighbors.
# Note for spatial models, n.samples is broken into a set of n.batch
# "batches", which each contain "batch.length" MCMC samples. In other words,
# n.samples = n.batch * batch.length
out.3 <- spPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2), 
                 det.formula = ~ scale(day) + (1 | obs),
                 data = data.veery, 
                 n.batch = 400, 
                 batch.length = 25,
		 priors = list(sigma.sq.unif = c(0, 20)),
		 NNGP = TRUE, 
		 n.neighbors = 10, 
                 n.thin = 10, 
                 n.burn = 5000, 
                 n.chains = 3,
                 n.report = 100)
summary(out.3)

# 3. Model validation -----------------------------------------------------
# Perform a posterior predictive check to assess model fit. 
ppc.out.1 <- ppcOcc(out.1, fit.stat = 'freeman-tukey', group = 1)
ppc.out.2 <- ppcOcc(out.2, fit.stat = 'freeman-tukey', group = 1)
ppc.out.3 <- ppcOcc(out.3, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out.1)
summary(ppc.out.2)
summary(ppc.out.3)

# 4. Model comparison -----------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
# Non-spatial
waicOcc(out.1)
# Non-spatial with intercept varying by BCR
waicOcc(out.2)
# Spatial model
waicOcc(out.3)

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out.3)
# Take a look at objects in resulting object
names(out.3)
str(out.3$beta.samples)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out.3$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out.3$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))

# 6. Prediction -----------------------------------------------------------
# Predict occupancy probability across the northeastern states
# Load prediction objects (loads objects elev.pred and coords.0)
# TODO: clearly something wrong with this prediction data here, as there is 
#       one value being repeated a billion times. 
load("data/bbs-2019-pred-data.rda")
# Standardize elevation prediction values by values used to fit model
elev.0 <- (elev.pred - mean(data.veery$occ.covs$elev)) / sd(data.veery$occ.covs$elev)
# Create prediction design matrix
X.0 <- cbind(1, elev.pred, elev.pred^2)
# Predict at new locations
out.pred <- predict(out.3, X.0)
# Occupancy probability means
psi.0.mean <- apply(out.pred$psi.0.samples, 2, mean)
psi.0.sd <- apply(out.pred$psi.0.samples, 2, sd)



# Create prediction map ---------------
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = my.proj)
# Full data
ne.states <- usa %>%
  dplyr::filter(ID %in% c('connecticut', 'delaware', 'maine', 'maryland',
		     'massachusetts', 'new hampshire', 'new jersey',
		     'new york', 'pennsylvania', 'rhode island',
		     'vermont'))
plot.df <- data.frame(psi.mean = psi.0.mean,
		      psi.sd = psi.0.sd,
		      x = coords.0[, 1],
		      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
ggplot(plot.df, aes(x = x, y = y, col = psi.mean)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_viridis_c()
psi.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.mean),interpolate = TRUE) +
  geom_sf(data = ne.states, alpha = 0) +
  scale_fill_gradientn("Richness Mean", colors = ocean.tempo(1000), limits = c(0, 1),
                       guide = guide_colourbar(title.position="top", reverse = FALSE)) +
  theme_bw(base_size = 25) +
  labs(x = "Longitude", y = "Latitude", fill = "") +
  theme(legend.position = c(0.85, 0.3),
        legend.background = element_rect(fill = NA))


