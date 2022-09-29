# dolphin-integrated-occupancy-model.R: this script fits an integrated occupancy model
#                                       using two independent data sources on bottlenose
#                                       dolphins in the northwest Mediterranean Sea from
#                                       Lauret et al. (2021) with the spOccupancy R package.
# Data source citation: 
#  Lauret, V., Labach, H., Authier, M., & Gimenez, O. (2021). 
#  Using single visits into integrated occupancy models to make the most of 
#  existing monitoring programs. Ecology, 102(12), e03535.

# Clear out the workspace.
rm(list = ls())
# Load the necessary packages.
library(spOccupancy)
# For plotting and summarizing results
library(MCMCvis)
library(ggplot2)
# If not using the RStudio project, set working directory to the repository
# directory. 
# setwd("../")

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.int)
load('data/lauret2021Ecology.rda')
# Check out the data list structure.  
str(data.int)
# Note the detection-nondetection data and the detection covariates are both 
# stored as lists, with separate elements for each data source. Within each list,
# the object is identical to that of single-species models. Additionally, 
# note the "site" element that defines which site the row of each data source
# corresponds to.


# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, integrated occupancy model. 
# Approximate run time: 1.2 minutes
out.full <- intPGOcc(occ.formula = ~ BATHY * SST,
		                 det.formula = list(samm = ~ eff.samm + ind.samm.2 + ind.samm.3 + ind.samm.4,
		         	                          gd = ~ eff.gd + ind.gd.2 + ind.gd.3 + ind.gd.4),
                     data = data.int,
                     n.samples = 10000, 
                     n.omp.threads = 1, 
                     verbose = TRUE, 
                     n.report = 2000, 
                     n.burn = 5000, 
                     n.thin = 10, 
                     n.chains = 3) 
summary(out.full)
# Fit non-spatial, integrated occupancy model with no interaction
# Approximate run time: 1.2 minutes
out.small <- intPGOcc(occ.formula = ~ BATHY + SST,
		                  det.formula = list(samm = ~ eff.samm + ind.samm.2 + ind.samm.3 + ind.samm.4,
		          	                         gd = ~ eff.gd + ind.gd.2 + ind.gd.3 + ind.gd.4),
                      data = data.int,
                      n.samples = 10000, 
                      n.omp.threads = 1, 
                      verbose = TRUE, 
                      n.report = 2000, 
                      n.burn = 5000, 
                      n.thin = 10, 
                      n.chains = 3) 

# 3. Model validation -----------------------------------------------------
# Perform a posterior predictive check to assess model fit. Note for integrated
# models a different posterior predictive check is given for each data set. 
ppc.out.full <- ppcOcc(out.full, fit.stat = 'freeman-tukey', group = 1)
ppc.out.small <- ppcOcc(out.small, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit.
summary(ppc.out.full)
summary(ppc.out.small)

# 4. Model comparison -----------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit. Note for integrated models a different
# WAIC value is given for each of the data sets.
waicOcc(out.full)
waicOcc(out.small)
# WAIC suggests the interaction is not that important (since the values are
# not that different between out.full and out.small), so we will continue
# working with the simpler model. 

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out.small)
# Take a look at objects in resulting object
names(out.small)
# Occupancy regression coefficient estimates
str(out.small$beta.samples)
# Probability that sea surface temperature (SST) has a positive value on occupancy
mean(out.small$beta.samples[, 3] > 0)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out.small$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out.small$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))

# 6. Prediction -----------------------------------------------------------
# Predict occupancy along a gradient of sea surface temperature.
# Create a set of values across the range of observed SST values (note
# values were standardized, so these are not on the real SST values cale).
SST.pred.vals <- seq(min(data.int$occ.covs$SST), max(data.int$occ.covs$SST), 
		                 length.out = 100)

# Predict occupancy across SST values at mean value of other occupancy covariates 
pred.df <- as.matrix(data.frame(intercept = 1, BATHY = 0, SST = SST.pred.vals))
out.pred <- predict(out.small, pred.df)
str(out.pred)
# Calculate the 2.5%, 50%, and 97.5% quantiles of the predicted occupancy values.
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile, prob = c(0.025, 0.5, 0.975))
# Format the data for plotting with ggplot2
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ],
                           psi.low = psi.0.quants[1, ],
                           psi.high = psi.0.quants[3, ],
                           SST = SST.pred.vals)
# Plot the predicted values
ggplot(psi.plot.dat, aes(x = SST, y = psi.med)) +
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), fill = 'grey70') +
  geom_line(size = 2) +
  theme_bw(base_size = 18) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Standardized SST', y = 'Occupancy Probability')

