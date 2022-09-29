# swiss-mhb-jsdm.R: this script fits a joint species distribution 
#                   model with imperfect detection via a factor
#                   modeling approach using data from the 
#                   Switzerland Breeding Bird Survey (Swiss MHB) 
#                   in 2014.
# Data source citations:   
# Kéry, M. & Royle, J.A. (2016) _Applied Hierarchical Modeling in Ecology_ AHM1 - 11.3.
# Swiss Federal Statistical Office (http://www.bfs.admin.ch)
# Data were derived from objects included in the AHMBook and unmarked R packages.
rm(list = ls())
library(spOccupancy)
# For summarizing results
library(MCMCvis)
# For making a species distribution map
library(ggplot2)
library(sf)
library(pals)
library(cowplot)
# For plotting residual covariance matrix
library(corrplot)
# If not using the RStudio project, set working directory to the repository
# directory. 
# setwd("../")
set.seed(101)

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.swiss.mhb)
load("data/swissMHB2014Data.rda")
# Check out the data list structure. Detection-nondetection data y
# are currently formatted for multiple species models. 
str(data.swiss.mhb)
# Subset the object to 20 species to speed things up
curr.sp <- sample(1:nrow(data.swiss.mhb$y), 20, replace = FALSE)
data.swiss.mhb$y <-data.swiss.mhb$y[curr.sp, , ]
str(data.swiss.mhb)

# 2. Model fitting --------------------------------------------------------
# Fit a spatially-explicit joint species distribution model with 
# imperfect detection. 
out <- sfMsPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest), 
                 det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur), 
                 data = data.swiss.mhb, 
                 n.batch = 400, 
                 batch.length = 25,
                 n.thin = 10, 
                 n.burn = 5000, 
                 n.chains = 1,
                 NNGP = TRUE,
                 n.factors = 4,
                 n.neighbors = 5,
                 n.omp.threads = 1,
                 cov.model = 'exponential',
                 n.report = 10)
summary(out, level = 'community')

# 3. Model validation -----------------------------------------------------
# Perform a posterior predictive check to assess model fit. 
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit. 
summary(ppc.out)

# 4. Model comparison -----------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. 
waicOcc(out)

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out, level = 'community')
summary(out, level = 'species')
summary(out, level = 'both')
# Create simple plot summaries using MCMCvis package.
# Occupancy community-level effects 
MCMCplot(out$beta.comm.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.comm.samples, ref_ovl = TRUE, ci = c(50, 95))
# Occupancy species-specific effects
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Species-specific factor loadings ----
# Useful to assess how many factor loadings to use in the model. If 
# the loadings (coefficients) for the last couple of factors are all very close
# to non-zero and significant, this suggests you can reduce the number of factors.
MCMCplot(out$lambda.samples, ref_ovl = TRUE, ci = c(50, 95))
# Species-species residual covariance matrix ---------
# Factor loadings 
lambda.means <- matrix(apply(out$lambda.samples, 2, mean), nrow(out$y), out$q)
lambda.means
# Calculate lambda %*% lambda as a singular species-species covariance matrix
Sigma.means <- lambda.means %*% t(lambda.means)
# Singular correlation matrix
species.cor.mat <- cov2cor(Sigma.means)
rownames(species.cor.mat) <- dimnames(data.swiss.mhb$y)[[1]]
colnames(species.cor.mat) <- dimnames(data.swiss.mhb$y)[[1]]
corrplot(species.cor.mat, method = 'square', type = 'lower')

# 6. Prediction -----------------------------------------------------------
# Predict occupancy probability and species richness across Switzerland
# Load prediction objects (loads objects pred.swiss and coords.0)
load("data/switzerlandPredData.rda")
str(pred.swiss)
# Note that prediction objects for multi-species models can be large, and 
# so instead of predicting at all 42,275 locations, we will predict at
# 8000 random non-sampled locations
set.seed(1000)
# Get the indices for the new locations (remove the locations where we have data)
new.site.indx <- which(is.na(match(do.call("paste", as.data.frame(coords.0)), 
                       do.call("paste", as.data.frame(data.swiss.mhb$coords)))))
pred.indx <- sample(new.site.indx, 8000, replace = FALSE)
pred.swiss <- pred.swiss[pred.indx, ]
coords.0 <- coords.0[pred.indx, ]
# plot(coords.0, pch = 19)
# Standardize elevation and forest prediction values by values used to fit model
elevation.0 <- (pred.swiss[, 'elevation'] - mean(data.swiss.mhb$occ.covs$elevation)) / 
                sd(data.swiss.mhb$occ.covs$elevation)
forest.0 <- (pred.swiss[, 'forest'] - mean(data.swiss.mhb$occ.covs$forest)) / 
                sd(data.swiss.mhb$occ.covs$forest)
# Create prediction design matrix
X.0 <- cbind(1, elevation.0, elevation.0^2, forest.0)
# Predict at new locations
out.pred <- predict(out, X.0, coords.0)
str(out.pred)

# Calculate species richness as a derived-quantity of the latent occupancy 
# values for each species
rich.samples <- apply(out.pred$z.0.samples, c(1, 3), sum)
# Mean species richness
rich.means <- apply(rich.samples, 2, mean)
# Standard deviation species richness
rich.sds <- apply(rich.samples, 2, sd)

# Create prediction map ---------------
plot.df <- data.frame(rich.mean = rich.means,
                      rich.sd = rich.sds,
                      x = coords.0[, 1], 
                      y = coords.0[, 2])
pred.sf <- st_as_sf(plot.df, coords = c('x', 'y'))
rich.mean.plot <- ggplot() +
  geom_sf(data = pred.sf, aes(color = rich.mean), size = 0.75) +
  scale_color_gradientn("Richness Mean", colors = ocean.tempo(1000)) +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "")

rich.sd.plot <- ggplot() +
  geom_sf(data = pred.sf, aes(color = rich.sd), size = 0.75) +
  scale_color_gradientn("Richness SD", colors = ocean.tempo(1000)) +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "")

plot_grid(rich.mean.plot, rich.sd.plot, ncol = 1)
