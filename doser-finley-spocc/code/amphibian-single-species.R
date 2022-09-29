# amphibian-single-species.R: this script fits a single-species occupancy model 
#                             using data from Ribeiro Jr. et al. (2018) with 
#                             the spOccupancy R package. 
# Data source citation:   
#   Ribeiro Jr, J. W., Siqueira, T., Brejão, G. L., & Zipkin, E. F. (2018). 
#   Effects of agriculture and topography on tropical amphibian species 
#   and communities. Ecological Applications, 28(6), 1554-1564.
rm(list = ls())
library(spOccupancy)
# For plotting and summarizing results
library(MCMCvis)
library(ggplot2)
# If not using the RStudio project, set working directory to the repository
# directory. 
# setwd("../")

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.crocar)
load("data/ribeiroJr2018EcoApps.rda")
str(data.crocar)

# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, single-species occupancy model. 
out <- PGOcc(occ.formula = ~ scale(forest) + scale(agriculture) + 
	                           scale(catchment) + scale(density) + 
	                           scale(slope), 
             det.formula = ~ scale(date) + I(scale(date)^2) + scale(rain), 
	           data = data.crocar, 
	           n.samples = 5000, 
	           n.thin = 2, 
	           n.burn = 3000, 
	           n.chains = 3,
	           n.report = 500)
summary(out)
# Fit a spatial, single-species occupancy model.
out.sp <- spPGOcc(occ.formula = ~ scale(forest) + scale(agriculture) + 
	                                scale(catchment) + scale(density) + 
	                                scale(slope), 
                  det.formula = ~ scale(date) + I(scale(date)^2) + scale(rain), 	          
	                data = data.crocar, 
		              n.batch = 400, 
		              batch.length = 25,
	                n.thin = 5, 
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
waicOcc(out)
waicOcc(out.sp)

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out)
# Take a look at objects in resulting object
names(out)
str(out$beta.samples)
# Probability the effect of stream density on occupancy is positive
mean(out$beta.samples[, 5] > 0)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects --------- 
MCMCplot(out$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))

# 6. Prediction -----------------------------------------------------------
# Predict occupancy along a gradient of forest cover.  
# Create a set of values across the range of observed forest values
forest.pred.vals <- seq(min(data.crocar$occ.covs$forest), 
			                  max(data.crocar$occ.covs$forest), 
			                  length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model
forest.pred.vals.scale <- (forest.pred.vals - mean(data.crocar$occ.covs$forest)) / 
	                         sd(data.crocar$occ.covs$forest)
# Predict occupancy across forest values at mean values of all other variables
pred.df <- as.matrix(data.frame(intercept = 1, forest = forest.pred.vals.scale, 
		                 agriculture = 0, catchment = 0, density = 0, 
		                 slope = 0))
out.pred <- predict(out, pred.df)
str(out.pred)
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile, 
		                  prob = c(0.025, 0.5, 0.975))
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ], 
			                     psi.low = psi.0.quants[1, ], 
			                     psi.high = psi.0.quants[3, ], 
                           forest = forest.pred.vals)
ggplot(psi.plot.dat, aes(x = forest, y = psi.med)) + 
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), fill = 'grey70') +
  geom_line() + 
  theme_bw() + 
  scale_y_continuous(limits = c(0, 1)) + 
  labs(x = 'Forest (% cover)', y = 'Occupancy Probability') 

