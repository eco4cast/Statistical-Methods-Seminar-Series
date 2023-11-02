# Date: 06 November 2023
# Author: Nicholas Clark (n.clark@uq.edu.au)

# Use Dynamic GAMs to analyze time series of desert rodent 
# relative abundances using long-term capture data from 
# the Portal Project

# Key references for the data:
# Ernest, SK Morgan, et al. 
#   "The Portal Project: a long-term study of a Chihuahuan 
#    desert ecosystem." BioRxiv (2018): 332783.

# Brown, James H., et al. 
#   "Complex species interactions and the dynamics 
#    of ecological systems: long-term experiments." 
#    Science 293.5530 (2001): 643-650.

####                  Load libraries                               ####
# This tutorial relies on the following packages:
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(dplyr)           # Tidy and flexible data manipulation
library(ggplot2)         # Flexible plotting
library(gratia)          # Graceful ggplot-based graphics for GAMs
library(marginaleffects) # Compute interpretable model predictions
library(tidybayes)       # Tidy manipulation / plots of posterior draws

# A custom ggplot2 theme; feel free to ignore if you have your own plot
# preferences
theme_set(theme_classic(base_size = 12, base_family = 'serif') +
            theme(axis.line.x.bottom = element_line(colour = "black",
                                                    size = 1),
                  axis.line.y.left = element_line(colour = "black",
                                                  size = 1)))
options(ggplot2.discrete.colour = c("#A25050",
                                    "#8F2727",
                                    'darkred',
                                    "#630000"),
        ggplot2.discrete.fill = c("#A25050",
                                  "#8F2727",
                                  'darkred',
                                  "#630000"))

####                  Inspect data                                 ####
# Load the pre-prepared Portal rodent abundance data
portal_ts <- read.csv('./data/portal_data.csv')

# Inspect the data structure, which contains lunar monthly total 
# captures across control plots for four different rodent species.
# It also contains a 12-month moving average of the unitless NDVI
# vegetation index, and monthly average minimum temperature
# (already scaled to unit variance)
dplyr::glimpse(portal_ts)

# The data contain 80 lunar monthly observations, though there are 
# plenty of NAs in the number of total captures (NAs are shown 
# as red bars in the below plot)
max(portal_ts$time)
image(is.na(t(portal_ts %>%
                dplyr::arrange(dplyr::desc(time)))), axes = F,
      col = c('grey80', 'darkred'))
axis(3, at = seq(0,1, len = NCOL(portal_ts)), 
     labels = colnames(portal_ts))

# The data are already in 'long' format, meaning each series x 
# time observation has its own entry in the dataframe.
# But {mvgam} needs a 'series' column that acts as a factor 
# indicator for the time series. Add one using {dplyr} commands:
portal_ts %>%
  dplyr::mutate(series = as.factor(species)) -> portal_ts
dplyr::glimpse(portal_ts)
levels(portal_ts$series)

# It is important that the number of levels matches the number of 
# unique series in the data to ensure indexing across series works 
# properly in the underlying modelling functions. For more information 
# on how to format data for {mvgam} modelling, see:
# https://nicholasjclark.github.io/mvgam/articles/data_in_mvgam.html

# Create a 'rel_abund' column that we can use as our outcome variable
# in a Beta regression model (see ?mvgam_families for more information
# on the kinds of observation models supported by {mvgam})
portal_ts %>%
  
  # Group by unique time points
  dplyr::group_by(time) %>%
  
  # Calculate total captures per time point
  dplyr::mutate(total = sum(captures, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  
  # Calculate relative abundance; {mvgam} does not yet allow zero- or 
  # one-inflated Beta observations, so we add small offsets for any
  # zeros and ones
  dplyr::mutate(rel_abund = 
                  pmin(0.999, pmax(0.001, 
                                   captures / total))) -> portal_ts
  
# Plot all of the time series together
plot_mvgam_series(data = portal_ts, y = 'rel_abund', series = 'all')

# Plot some more in-depth features for individual series
plot_mvgam_series(data = portal_ts, y = 'rel_abund', series = 1)
plot_mvgam_series(data = portal_ts, y = 'rel_abund', series = 2)
plot_mvgam_series(data = portal_ts, y = 'rel_abund', series = 3)
plot_mvgam_series(data = portal_ts, y = 'rel_abund', series = 4)

# Inspect some associations between logit(rel_abund) and 
# minimum temperature / NDVI moving average for each species 
# to get a sense of how their relative abundances vary over seasons 
# and with varying habitat conditions

# DM
portal_ts %>% 
  dplyr::filter(species == 'DM') %>%
  ggplot(aes(x = mintemp, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'DM',
       y = "logit(relative abundance)", 
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'DM') %>%
  ggplot(aes(x = ndvi_ma12, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')

# PP
portal_ts %>% 
  dplyr::filter(species == 'PP') %>%
  ggplot(aes(x = mintemp, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'PP',
       y = "logit(relative abundance)", 
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'PP') %>%
  ggplot(aes(x = ndvi_ma12, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')

# PB
portal_ts %>% 
  dplyr::filter(species == 'PB') %>%
  ggplot(aes(x = mintemp, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'PB',
       y = "logit(relative abundance)",
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'PB') %>%
  ggplot(aes(x = ndvi_ma12, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')

# DO
portal_ts %>% 
  dplyr::filter(species == 'DO') %>%
  ggplot(aes(x = mintemp, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(title = 'DO',
       y = "logit(relative abundance)", 
       x = 'Minimum temperature') +
  
  portal_ts %>% 
  dplyr::filter(species == 'DO') %>%
  ggplot(aes(x = ndvi_ma12, y = qlogis(rel_abund))) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 10),
              col = 'darkred', fill = "#A25050") +
  labs(y = NULL, 
       x = 'NDVI moving average')

# There may be some support for some nonlinear effects here; 
# Now on to some modelling. But first we will split the data 
# into training and testing folds to evaluate predictions
data_train <- portal_ts %>%
  dplyr::filter(time <= 68)
data_test <- portal_ts %>%
  dplyr::filter(time > 68)

####                  Fit DGAMs                                   ####
# An initial model will attempt to capture variation 
# in species' responses to minimum temperature, using a 
# Hierarchical GAM (HGAM) with no dynamic component 
# (see ?mgcv::s, ?mgcv::te, ?brms::gp and ?mvgam::dynamic for 
# more information on the kinds of terms supported in {mvgam}
# formulae). This model takes ~ 20 seconds to fit, after compilation
mod <- mvgam(rel_abund ~ 
               
               # Hierarchical intercepts capture variation in average
               # relative abundances
               s(series, bs = 're') +
               
               # A shared smooth of minimum temperature
               s(mintemp, k = 8) +
               
               # Deviation smooths of minimum temperature,
               # allowing each species' response to mintemp to vary
               # from the shared smooth
               s(mintemp, series, bs = 'sz', k = 8) - 1,
             
             # Condition on the training data
             data = data_train,
             
             # Automatically compute forecasts for the test data
             newdata = data_test,
             
             # Beta observations with independent precisions
             family = betar(),
             
             # cmdstanr is highly recommended over rstan as 
             # it is much more up-to-date with the latest 
             # features in the Stan universe
             backend = 'cmdstanr')

# The model can be described mathematically as follows:

#                   for s in 1:N_species ...
#                  for t in 1:N_timepoints...

#                   ## Observation model ##
#  Relative abundance[s, t] ~ Beta(μ[s, t], φ[s])

#                   ## Linear predictor ##
#            logit(μ[s, t]) = α[s] + f(mintemp)_shared[t] + 
#                             f(mintemp)_species[s, t]
#                         f = sum(β_smooth * b) 

#                      ## Priors ##
#                         α ~ Normal(μ_population, σ_population)
#              μ_population ~ Normal(0, 1)
#              σ_population ~ Student-t(3, 0, 2.5)[0, ]
#                  β_smooth ~ MVNormal(0, (Ω ∗ λ)^−1)
#                         λ ~ Normal(5, 30)[0, ]
#                         φ ~ Gamma(0.01, 0.01)

# where: 
# f are the penalized smooth functions
# b are thin plate spline basis functions
# Ω are a set of prior precision matrices for the smooths
# λ are the smoothing penalties that prevent overfitting; note that
#   Normal(5, 30)[0, ] indicates a half-normal prior distribution
# φ are species-level precision parameters

# If you would like to see the underlying Stan code, which is fully
# transparent in its use of prior distributions, use the code()
# function:
code(mod)

# Inspect the model summary
summary(mod)

# Sampling diagnostics (see ?mcmc_plot.mvgam for details on the types
# of {bayesplot} plots that can be used with {mvgam})
mcmc_plot(mod, type = 'rhat_hist')
mcmc_plot(mod, type = 'trace')

# Pairs plots are also useful for diagnosing non-identifiabilities in
# Bayesian models. See ?bayesplot::mcmc_pairs for details. Here a pairs plot
# of the random effect mean and SD parameters shows no worrisome 'funnel' 
# behaviour that can plague hierarchical models:
pairs(mod, variable = c('mean(series)', 'sd(series)'))

# Plot the hierarchical smooth components with the S3 'plot' function
# (see ?plot.mvgam for more details)
plot(mod, type = 'smooths')

# Plot the hierarchical intercepts
plot(mod, type = 're')

# More informative plots can be made using plot_predictions() from
# the {marginaleffects} universe to visualise conditional effects 
# on the outcome scale. See ?marginaleffects::plot_predictions 
# for details, or visit: https://marginaleffects.com/
plot_predictions(mod, 
                 condition = c('mintemp', 'series', 'series'),
                 points = 0.5,
                 rug = TRUE) +
  theme(legend.position = 'none') +
  labs(y = 'Relative abundance', x = 'Minimum temperature')

# A similar plot can be made on the link scale, which shows that 
# the temperature effect is mostly estimated to be linear
plot_predictions(mod, 
                 condition = c('mintemp', 'series', 'series'),
                 type = 'link') +
  theme(legend.position = 'none') +
  labs(y = 'logit(relative abundance)', x = 'Minimum temperature')

# The average (marginal) effect of mintemp can be plotted for each
# series using the plot_slopes() function:
plot_slopes(mod, variable = 'mintemp',
            condition = c('series', 'series'),
            type = 'link') +
  theme(legend.position = 'none') +
  labs(y = 'logit(relative abundance)', x = 'Series')

# Predictions are your best tool when trying to interpret / report
# GAMs (or any regression-based model for that matter). As an 
# illustration, we can compute how much higher we would expect 
# relative abundance to be if the scaled mintemp variable was 
# high (1.75) vs low (-1.75) for each species. 
# We can take full advantage of the comparisons function 
# from the {marginaleffects} package to do this 
# (see ?marginaleffects::comparisons for details)
post_contrasts <- comparisons(mod,
                              
                              # Set the contrast to compute
                              variables = 
                                list(mintemp = c(-1.75, 1.75)),
                              
                              # Compute contrast for each level of
                              # the 'series' variable
                              newdata = 
                                datagrid(series = unique)) %>%
  
  # Take posterior draws
  posteriordraws()

# Use the resulting posterior draw object to plot a density of the 
# posterior contrasts
post_contrasts %>%
  ggplot(aes(x = draw)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  
  # Use the stat_halfeye function from {tidybayes} for a nice visual
  stat_halfeye(fill = "#C79999", alpha = 0.75) +
  facet_wrap(~ series, scales = 'free') +
  labs(x = "Change in relative abundance if mintemp is high vs low",
       y = "Density", 
       title = "Posterior temperature contrasts") +
  theme_classic()

# How do the posterior predictions look for the hindcast 
# and forecast periods?
plot(mod, type = 'forecast', series = 1)
plot(mod, type = 'forecast', series = 2)
plot(mod, type = 'forecast', series = 3)
plot(mod, type = 'forecast', series = 4)

# Randomized Quantile (Dunn-Smyth) residuals provide model diagnostics
# Dunn, Peter K., and Gordon K. Smyth. 
#    "Randomized quantile residuals." 
#    Journal of Computational and graphical statistics 5.3 
#    (1996): 236-244.
plot(mod, type = 'residuals', series = 1)
plot(mod, type = 'residuals', series = 2)
plot(mod, type = 'residuals', series = 3)
plot(mod, type = 'residuals', series = 4)

# We are clearly missing a lot of the variation in the data. 
# It makes sense to think about interactions between temperature
# and NDVI in this system, as the green-up period happens at 
# certain times of year. We can incorporate this knowledge 
# into a more complex HGAM (takes ~ 60 seconds to fit)
mod2 <- mvgam(rel_abund ~ 
               
               # Hierarchical intercepts capture variation in average
               # relative abundances
               s(series, bs = 're') +
                
               # Hierarchical tensor products of mintemp and 
               # NDVI capture nonlinear interactions of these two
               # covariates
               te(mintemp, ndvi_ma12, k = c(3, 5)) +
               te(mintemp, ndvi_ma12, by = series,
                  k = c(3, 5), m = 1) - 1,
             
             data = data_train,
             newdata = data_test,
             family = betar(),
             backend = 'cmdstanr')

# The model is not much different to the one above, and 
# can be described as follows:

#                   for s in 1:N_species ...
#                  for t in 1:N_timepoints...

#                   ## Observation model ##
#  Relative abundance[s, t] ~ Beta(μ[s, t], φ[s])

#                   ## Linear predictor ##
#            logit(μ[s, t]) = α[s] + f(mintemp, ndvi)_shared[t] + 
#                             f(mintemp, ndvi)_species[s, t]
#                         f = sum(β_smooth * b) 

#                      ## Priors ##
#                         α ~ Normal(μ_population, σ_population)
#              μ_population ~ Normal(0, 1)
#              σ_population ~ Student-t(3, 0, 2.5)[0, ]
#                  β_smooth ~ MVNormal(0, (Ω ∗ λ)^−1)
#                         λ ~ Normal(5, 30)[0, ]
#                         φ ~ Gamma(0.01, 0.01)

# The model summary now has many effects in it; 
# suppress the un-interpretable spline coefficients
summary(mod2, include_betas = FALSE)

# Sampling diagnostics
mcmc_plot(mod2, type = 'rhat_hist')

# The hierarchical smooth components are now a series of bivariate
# nonlinear interaction smooths
plot(mod2, type = 'smooths')

# {gratia} makes these a bit easier to visualise
gratia::draw(mod2$mgcv_model)

# But as before, we can't get a good handle on what 
# these effects mean without some useful conditional predictions;
# I've reduced the confidence level to 50% here so the uncertainty
# intervals are a bit easier to see
plot_predictions(mod2, 
                 condition = c('mintemp', 'ndvi_ma12', 'series'),
                 points = 0.5, conf_level = 0.5) +
  labs(y = 'Relative Abundance', x = 'Minimum temperature')

plot_predictions(mod2, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'),
                 points = 0.5, conf_level = 0.5) +
  labs(y = 'Relative abundance', x = 'NDVI moving average')

# As before we can also plot these on the link scale
plot_predictions(mod2, 
                 condition = c('mintemp', 'ndvi_ma12', 'series'),
                 type = 'link') +
  labs(y = 'logit(relative abundance)', x = 'Minimum temperature')

plot_predictions(mod2, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'),
                 type = 'link') +
  labs(y = 'logit(relative abundance)', x = 'NDVI moving average')

# In-sample fit metrics using approximate leave-one-out cross-validation
# slightly prefer the more complex model, but there is large 
# uncertainty (see ?loo::loo for details on  how PSIS-LOO is calculated)
loo(mod)$estimates
loo(mod2)$estimates
loo_compare(mod, mod2)

# Posterior predictive checks can be another useful way to ask
# whether the model is able to reproduce key features of the observed
# data. Here we plot predictive histograms and Probability Integral
# Transform (PIT) histograms to show that both models produce similar
# predictions for the first rodent species (DM). See ?ppc for more details
# on the kinds of predictive checks supported by mvgam
layout(matrix(1:4, nrow = 2, byrow = TRUE))
ppc(mod, type = 'hist', series = 1)
title(main = 'Model 1')
ppc(mod, type = 'pit', series = 1)

ppc(mod2, type = 'hist', series = 1)
title(main = 'Model 2')
ppc(mod2, type = 'pit', series = 1)
layout(1)

# Posterior predictions still leave a lot to be desired
plot(mod2, type = 'forecast', series = 1)
plot(mod2, type = 'forecast', series = 2)
plot(mod2, type = 'forecast', series = 3)
plot(mod2, type = 'forecast', series = 4)

# And Randomized Quantile residuals show clear patterns 
# of unmodelled temporal autocorrelation
plot(mod2, type = 'residuals', series = 1)
plot(mod2, type = 'residuals', series = 2)
plot(mod2, type = 'residuals', series = 3)
plot(mod2, type = 'residuals', series = 4)

# We need to capture this autocorrelation, which is best done using 
# a latent dynamic process. Here we build a State-Space model using
# the 'trend_formula' argument. This model assumes that environmental
# predictors act on the 'true' but latent relative abundance for each
# species, but there is also strong autocorrelation that needs to
# be captured (which we model with an AR1 process). The observation
# formula is empty, but we could incorporate other predictors here
# if we thought they might influence our ability to take accurate
# measurements (such as storms, unusually cold nights, trap failures
# etc...; takes ~ 1 minute to fit after compilation)
mod3 <- mvgam(rel_abund ~ -1,
                
              # GAM components now moved to the latent process model
              trend_formula = ~
                
                # Hierarchical intercepts capture variation in average
                # relative abundances (we use 'trend' rather
                # than 'series' here; this is because it is possible 
                # to set up models in which some series share the same
                # latent dynamic process, so we can have fewer process
                # models than we have series). See examples of this here:
                # https://nicholasjclark.github.io/mvgam/articles/shared_states.html
                s(trend, bs = 're') +
              
                # Hierarchical tensor products of mintemp and 
                # NDVI capture nonlinear interactions of these two
                # covariates
                te(mintemp, ndvi_ma12, k = c(3, 5)) +
                te(mintemp, ndvi_ma12, by = trend,
                   k = c(3, 5), m = 1),
              
              # Independent latent AR1 processes; see ?mvgam_trends
              # for details on the types of dynamic processes currently
              # supported by {mvgam}
              trend_model = 'AR1',
              
              # Placing more informative priors on the AR1 
              # and process error parameters can be useful to incorporate
              # our prior domain expertise. see ?get_mvgam_priors
              # and ?brms::prior for more information about how to 
              # see what parameters can have priors changed and how 
              # to change priors
              priors = c(prior(exponential(1), class = sigma),
                         prior(normal(0.5, 0.25), class = ar1,
                             lb = -1, ub = 1)),
              
              data = data_train,
              newdata = data_test,
              family = betar(),
              backend = 'cmdstanr')

# The model introduces latent process models, and 
# can be described as follows:

#                   for s in 1:N_species ...
#                  for t in 1:N_timepoints...

#                   ## Observation model ##
#  Relative abundance[s, t] ~ Beta(μ[s, t], φ[s])

#               ## Observation linear predictor ##
#            logit(μ[s, t]) = z[s, t]

#               ## Dynamic process model ##
#                   z[s, t] ~ Normal(μ_process[s, t], σ_process[s])
#           μ_process[s, t] = α[s] + f(mintemp, ndvi)_shared[t] + 
#                             f(mintemp, ndvi)_species[s, t] + 
#                             AR1 * μ_process[s, t - 1]
#                         f = sum(β_smooth * b) 

#                      ## Priors ##
#                         α ~ Normal(μ_population, σ_population)
#              μ_population ~ Normal(0, 1)
#              σ_population ~ Student-t(3, 0, 2.5)[0, ]
#                  β_smooth ~ MVNormal(0, (Ω ∗ λ)^−1)
#                         λ ~ Normal(5, 30)[0, ]
#                       AR1 ~ Normal(0.5, 0.25)[-1, 1]
#                 σ_process ~ Exponential(1)
#                         φ ~ Gamma(0.01, 0.01)
# where: 
# σ_process is the process error for each species
# AR1 captures temporal autocorrelation in the process models

# Again you can inspect the Stan code to see how this is all set up
code(mod3)

# The model summary now includes summaries of the latent 
# process parameters
summary(mod3, include_betas = FALSE)

# The AR1 parameters are all strongly positive, suggesting they play
# a key role in capturing unmodelled temporal dynamics
mcmc_plot(mod3, variable = 'ar1', regex = TRUE, type = 'areas')

# Sampling diagnostics are still reassuring
mcmc_plot(mod3, type = 'rhat_hist')

# The hierarchical smooth components can now be visualised using
# 'trend_effects = TRUE'
plot(mod3, type = 'smooths', trend_effects = TRUE)
layout(1)

# plot_predictions() works as before
plot_predictions(mod3, 
                 condition = c('mintemp', 'ndvi_ma12', 'series'),
                 points = 0.5) +
  labs(y = 'Relative abundance', x = 'Minimum temperature')

plot_predictions(mod3, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'),
                 points = 0.5) +
  labs(y = 'Relative abundance', x = 'NDVI moving average')

# And again on the link scale
plot_predictions(mod3, 
                 condition = c('mintemp', 'ndvi_ma12', 'series'),
                 type = 'link') +
  labs(y = 'logitaRelative abundance)', x = 'Minimum temperature')

plot_predictions(mod3, 
                 condition = c('ndvi_ma12', 'mintemp', 'series'),
                 type = 'link') +
  labs(y = 'logit(relative abundance)', x = 'NDVI moving average')

# Hindcasts and forecasts are much better, as we would expect
plot(mod3, type = 'forecast', series = 1)
plot(mod3, type = 'forecast', series = 2)
plot(mod3, type = 'forecast', series = 3)
plot(mod3, type = 'forecast', series = 4)

# Residual diagnostics now look better as well
plot(mod3, type = 'residuals', series = 1)
plot(mod3, type = 'residuals', series = 2)
plot(mod3, type = 'residuals', series = 3)
plot(mod3, type = 'residuals', series = 4)

# Of course this model has better in-sample predictive performance
loo_compare(mod3, mod2, mod)

# PPCs of predictive densities also look reasonable
layout(matrix(1:4, nrow = 2, byrow = TRUE))
ppc(mod3, type = 'density', series = 1)
ppc(mod3, type = 'density', series = 2)
ppc(mod3, type = 'density', series = 3)
ppc(mod3, type = 'density', series = 4)
layout(1)

# Evaluating forecasts using proper scoring rules is an important step
# to compare different models. {mvgam} can score objects of class
# mvgam_forecast. See ?forecast.mvgam and ?score.mvgam_forecast for
# more details
fc <- forecast(mod3, newdata = data_test)
class(fc)

# Compute the energy score, which targets calibration of multivariate
# forecast distributions against out of sample observations
energy_fc <- score(fc, score = 'energy')
str(energy_fc)
energy_fc$DM
energy_fc$all_series

# Comparing these scores among different candidate models is helpful
# to decide on which model(s) to use in further applications
# Plot the difference in energy scores for this model and the above
# (mod2); a negative value means the State-Space model is better,
# while a positive value means the non-dynamic model is better
par(mar = c(3.5, 3.5, 2, 2),
    mgp = c(2, 0.5, 0))
diff_scores <- energy_fc$all_series$score -
  score(forecast(mod2), score = 'energy')$all_series$score
plot(diff_scores, pch = 16, col = 'darkred', 
     ylim = c(-1*max(abs(diff_scores), na.rm = TRUE),
              max(abs(diff_scores), na.rm = TRUE)),
     bty = 'l',
     xlab = 'Forecast horizon',
     ylab = expression(Energy[dynamic]~-~Energy[non-dynamic]))
box(bty = 'l', lwd = 2)
abline(h = 0, lty = 'dashed', lwd = 2)

####                                                                 ####
####                      Not run in webinar                         ####
####                                                                 ####

# Below are some examples to show how we can incorporate 
# multivariate dynamics. These won't be run in the webinar 
# as they take a bit longer to fit, but they can be useful examples 
# if you are interested in capturing dependencies in collections 
# of ecological time series. 

# See more about multivariate dynamics in {mvgam} here: 
# https://nicholasjclark.github.io/mvgam/articles/trend_formulas.html
# and here: 
# https://nicholasjclark.github.io/physalia-forecasting-course/day4/tutorial_4_physalia

# Dynamic factor processes, where a set of dimension-reduced
# factors is propagated over time, and each series is allowed to depend
# on these factors via a set of loading coefficients that must be
# estimated
# Takes ~ 200 seconds to fit, after compilation
mod4 <- mvgam(rel_abund ~ 
                
                # Hierarchical intercepts capture variation in average
                # relative abundances
                s(series, bs = 're') +
                
                # Hierarchical tensor products of mintemp and 
                # NDVI capture nonlinear interactions of these two
                # covariates
                te(mintemp, ndvi_ma12, k = c(3, 5)) +
                te(mintemp, ndvi_ma12, by = series,
                   k = c(3, 5), m = 1) - 1,
              
              # A set of dynamic factors evolve as Gaussian Processes
              # to capture the temporal dependencies in the data
              trend_model = 'GP',
              use_lv = TRUE,
              n_lv = 2,
              
              # Usual data and response distribution
              data = data_train,
              newdata = data_test,
              family = betar(),
              backend = 'cmdstanr')

# The usual summaries
summary(mod4, include_betas = FALSE)
mcmc_plot(mod4, type = 'rhat_hist')

# Plot the dynamic factors
plot_mvgam_factors(mod4)

# Loadings on these factors can induce dependence among the
# time series; visualize these dependencies as correlations
correlations <- lv_correlations(object = mod4)
mean_correlations <- correlations$mean_correlations
mean_correlations[upper.tri(mean_correlations)] <- NA
mean_correlations <- data.frame(mean_correlations)
mean_correlations %>%
  dplyr::add_rownames("series1") %>%
  tidyr::pivot_longer(-c(series1), 
                      names_to = "series2", 
                      values_to = "Correlation") -> mean_correlations
ggplot(mean_correlations,
       aes(x = series1, y = series2)) + 
  geom_tile(aes(fill = Correlation)) +
  scale_fill_gradient2(low = "darkred", 
                       mid = "white", 
                       high = "darkblue",
                       midpoint = 0,
                       breaks = seq(-1,1,length.out = 5),
                       limits = c(-1, 1),
                       name = 'Trend\ncorrelation') + 
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# Notice how some of the series have very similar (correlated)
# dynamic trend components
plot(mod4, type = 'trend', series = 1)
plot(mod4, type = 'trend', series = 2)
plot(mod4, type = 'trend', series = 3)
plot(mod4, type = 'trend', series = 4)

# This also forces them to have correlated forecasts
plot(mod4, type = 'forecast', series = 1)
plot(mod4, type = 'forecast', series = 2)
plot(mod4, type = 'forecast', series = 3)
plot(mod4, type = 'forecast', series = 4)

# Multivariate autoregressive processes
# This model fits an equivalent model to mod3 above, except it uses
# a VAR1 for the dynamics rather than independent AR1 processes. This 
# can capture lagged cross-dependence among the process models, and can 
# also capture correlated process errors
# Takes ~ 8 minutes to fit, after compilation
mod5 <- mvgam(rel_abund ~ -1,
              
              # GAM components now moved to the latent process model
              trend_formula = ~
                
                # Hierarchical intercepts capture variation in average
                # relative abundances (we use 'trend' rather
                # than 'series' here; this is because it is possible 
                # to set up models in which some series share the same
                # latent dynamic process, so we can have fewer process
                # models than we have series)
                s(trend, bs = 're') +
                
                # Hierarchical tensor products of mintemp and 
                # NDVI capture nonlinear interactions of these two
                # covariates
                te(mintemp, ndvi_ma12, k = c(3, 5)) +
                te(mintemp, ndvi_ma12, by = trend,
                   k = c(3, 5), m = 1),
              
              # Vector Autoregressive dynamics capture lagged
              # cross-dependence and correlated process errors
              trend_model = 'VAR1cor',
              
              # Prior for process error variance components 
              # (the covariances are parameterized by a default
              # LKJ prior on the error correlations)
              priors = prior(exponential(1), class = sigma),
              
              # Usual data and response arguments
              data = data_train,
              newdata = data_test,
              family = betar(),
              backend = 'cmdstanr')

# This summary is more complex than the others, as we get posterior
# estimates for the full VAR1 coefficient matrix (labelled as 'A') and 
# the full process error covariance matrix ('Sigma')
summary(mod5, include_betas = FALSE)

# View estimates for the VAR1 matrix to see there
# is some support, but not much for lagged cross-dependence
# on the off-diagonals. Each cell captures the lagged effect 
# of the process in the column on the process in the row for 
# the next timestep (i.e. cell [3,1] captures effect of species 1's
# process estimate on the process estimate for species 3 in the 
# next timestep)
A_pars <- matrix(NA, nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    A_pars[i, j] <- paste0('A[', i, ',', j, ']')
  }
}
mcmc_plot(mod5, 
          variable = as.vector(t(A_pars)), 
          type = 'hist')


# View estimates for the process error covariance matrix to see there
# is support for some correlated process errors on the off-diagonals
Sigma_pars <- matrix(NA, nrow = 5, ncol = 5)
for(i in 1:5){
  for(j in 1:5){
    Sigma_pars[i, j] <- paste0('Sigma[', i, ',', j, ']')
  }
}
mcmc_plot(mod5, 
          variable = as.vector(t(Sigma_pars)), 
          type = 'hist')

# Comparing energy scores for this model and the
# AR1 model (without multispecies dependence) suggests the 
# simpler model does better in some of the validation timepoints, and
# the differences in scores tend to be very small.
# This is ikely because there was a considerable regime 
# shift happening during the forecast horizon, leading to different 
# dynamics compared to what was observed in the training data and 
# perhaps favouring simpler models
diff_scores <- score(forecast(mod5), 
                     score = 'energy')$all_series$score -
  score(forecast(mod3), 
        score = 'energy')$all_series$score
plot(diff_scores, pch = 16, col = 'darkred', 
     ylim = c(-1*max(abs(diff_scores), na.rm = TRUE),
              max(abs(diff_scores), na.rm = TRUE)),
     bty = 'l',
     xlab = 'Forecast horizon',
     ylab = expression(Energy[VAR1]~-~Energy[AR1]))
box(bty = 'l', lwd = 2)
abline(h = 0, lty = 'dashed', lwd = 2)
