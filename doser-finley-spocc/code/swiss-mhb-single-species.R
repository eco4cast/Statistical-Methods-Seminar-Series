# swiss-mhb-single-species.R: this script fits a single-species occupancy model 
#                             using data on the European Goldfinch
#                             from the Switzerland Breeding Bird Survey 
#                             (Swiss MHB) in 2014. 
# Data source citations:   
# Kéry, M. & Royle, J.A. (2016) _Applied Hierarchical Modeling in Ecology_ AHM1 - 11.3.
# Swiss Federal Statistical Office (http://www.bfs.admin.ch)
# Data were derived from objects included in the AHMBook and unmarked R packages.
rm(list = ls())
library(spOccupancy)
# For summarizing MCMC results
library(MCMCvis)
# For making species distribution maps
library(ggplot2)
library(stars)
library(pals)
library(cowplot)
# If not using the RStudio project, set working directory to the repository
# directory. 
# setwd("../")
# Set seed for same results
set.seed(250)

# In this example, our goal is to produce a species distribution map of the 
# European Goldfinch throughout Switzerland.

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.goldfinch)
load("data/europeanGoldfinchSwiss.rda")
str(data.goldfinch)
# Take a quick look at the spatial locations
plot(data.goldfinch$coords, pch = 19)

# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, single-species occupancy model
out <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest), 
             det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur), 
             data = data.goldfinch, 
             n.samples = 5000, 
             n.thin = 4, 
             n.burn = 3000, 
             n.chains = 3,
             n.report = 500)
summary(out)

# Fit a spatial, single-species occupancy model using an NNGP and 10 neighbors
# Note for spatial models, n.samples is broken into a set of "n.batch"
# batches, which each contain "batch.length" MCMC samples. In other words,
# n.samples = n.batch * batch.length
out.sp <- spPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest), 
                  det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur), 
                  data = data.goldfinch, 
                  n.batch = 400, 
                  batch.length = 25,
                  NNGP = TRUE, 
                  n.neighbors = 5, 
                  n.thin = 10, 
                  n.burn = 5000, 
                  n.chains = 3,
                  n.report = 100)
summary(out.sp)

# 3. Model validation -----------------------------------------------------
# Perform a posterior predictive check to assess model fit. 
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
ppc.out.sp <- ppcOcc(out.sp, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out)
summary(ppc.out.sp)

# 4. Model comparison -----------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
# Non-spatial
waicOcc(out)
# Spatial
waicOcc(out.sp)

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out.sp)
# Take a look at objects in resulting object
names(out.sp)
str(out.sp$beta.samples)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out.sp$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out.sp$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))

# 6. Prediction -----------------------------------------------------------
# Predict occupancy probability across Switzerland
# Load prediction objects (loads objects pred.swiss and coords.0)
load("data/switzerlandPredData.rda")
str(pred.swiss)
# Standardize elevation and forest prediction values by values used to fit model
elevation.0 <- (pred.swiss[, 'elevation'] - mean(data.goldfinch$occ.covs$elevation)) / 
                sd(data.goldfinch$occ.covs$elevation)
forest.0 <- (pred.swiss[, 'forest'] - mean(data.goldfinch$occ.covs$forest)) / 
                sd(data.goldfinch$occ.covs$forest)
# Create prediction design matrix
X.0 <- cbind(1, elevation.0, elevation.0^2, forest.0)
# Predict at new locations
out.pred <- predict(out.sp, X.0, coords.0)
# Occupancy probability means
psi.0.mean <- apply(out.pred$psi.0.samples, 2, mean)
# Occupancy probability standard deviations
psi.0.sd <- apply(out.pred$psi.0.samples, 2, sd)
# Spatial process mean and sd
w.0.mean <- apply(out.pred$w.0.samples, 2, mean)
w.0.sd <- apply(out.pred$w.0.samples, 2, sd)


# Create a species distribution map with uncertainty ----------------------
plot.df <- data.frame(psi.mean = psi.0.mean,
                      psi.sd = psi.0.sd,
                      w.mean = w.0.mean, 
                      w.sd = w.0.sd,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
psi.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.mean),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy Mean')
psi.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.sd),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy SD')
w.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.mean),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Spatial Effect Mean')
w.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.sd),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Spatial Effect SD') 
plot_grid(psi.mean.plot, w.mean.plot, 
          psi.sd.plot, w.sd.plot, nrow = 2, ncol = 2)
