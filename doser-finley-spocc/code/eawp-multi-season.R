# eawp-multi-season.R: this script fits a multi-season occupancy model using
#                      a small subset of data from Doser et al. (2021) 
#                      using the spOccupancy R package.
# Data source citation: 
#   Doser, J. W., Weed, A. S., Zipkin, E. F., Miller, K. M., & Finley, A. O. (2021). 
#     Trends in bird abundance differ among protected forests but not bird guilds. 
#     Ecological Applications, 31(6), e02377.
#   Faccio, S., B. R. Mitchell, and P. S. Pooler. 2015. 
#     Northeast Temperate Network breeding landbird monitoring protocol: 2015 revision. 
#     Technical Report Natural Resource Report NPS/NETN/NRR-2015/942. 
#     National Park Service, Fort Collins, Colorado, USA.

rm(list = ls())
library(spOccupancy)
# For plotting and summarizing results
library(MCMCvis)
library(ggplot2)
# If not using the RStudio project, set working directory to the repository
# directory. 
# setwd("../")

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.list)
load("data/doser2021EcoApps.rda")
str(data.list)
# Note the detection-nondetection data are stored as a three-dimensional
# array with dimensions of site, year, and replication. Also note the 
# occupancy covariates are now stored as a list as they can vary across both 
# space and season.

# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, multi-season occupancy model. 
out <- tPGOcc(occ.formula = ~ regen + basalArea + I(basalArea^2) + percentForest + 
                              (1 | site.index) + scale(year),
              det.formula = ~ basalArea + (1 | year), 
	            data = data.list, 
	            n.batch = 200,
	            batch.length = 25,
	            ar1 = TRUE,
	            n.thin = 2, 
	            n.burn = 3000, 
	            n.chains = 3,
	            n.report = 50)
summary(out)
# Fit a spatial, multi-season occupancy model.
out.sp <- stPGOcc(occ.formula = ~ regen + basalArea + I(basalArea^2) + percentForest + (1 | site.index) + 
	                                scale(year),
                  det.formula = ~ basalArea + (1 | year), 
	                data = data.list, 
	                n.batch = 200,
	                batch.length = 25,
	                ar1 = TRUE,
	                n.thin = 2, 
	                n.burn = 3000, 
	                n.chains = 3,
	                n.report = 50)
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
waicOcc(out)
waicOcc(out.sp)

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out)
# Take a look at objects in resulting object
names(out)
str(out$beta.samples)
# Probability the linear temporal trend is negative
mean(out$beta.samples[, 6] < 0)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
# Compuate average occupancy probability across all 25 sites during each year
psi.avg.by.year <- apply(out$psi.samples, 3, mean)
# Compute 95% credible interval of average values during each year
psi.ci.by.year <- apply(out$psi.samples, 3, quantile, c(0.025, 0.975))
plot.df <- data.frame(year = min(data.list$occ.covs$year):max(data.list$occ.covs$year), 
                      psi.mean = psi.avg.by.year, 
                      psi.low = psi.ci.by.year[1, ], 
                      psi.high = psi.ci.by.year[2, ])
ggplot(data = plot.df, aes(x = year, y = psi.mean)) + 
  geom_point(size = 2.5) + 
  geom_segment(aes(x = year, y = psi.low, yend = psi.high, xend = year), size = 1) + 
  theme_bw(base_size = 18) + 
  labs(x = 'Year', y = 'Average Occupancy Probability')
