#####-----Advances in piecewise estimation of path models-----------------------

# EFI/ESA Seminar
# Date: 05 December 2022
# Author: Jon Lefcheck
# Contact: LefcheckJ@si.edu

# Load required libraries
# devtools::install_github("jslefche/piecewiseSEM@devel") # version 2.3.0
library(piecewiseSEM)
library(lavaan)

# Load Keeley data set
data(keeley)

# Examine Keeley data
head(keeley)

#####-----Fit structural equation model to Keeley data--------------------------

# Fit simple/multiple regressions using `lm`

# Break down component regressions
abiotic_model <- lm(abiotic ~ distance, data = keeley)
hetero_model <- lm(hetero ~ distance, data = keeley)
richness_model <- lm(rich ~ abiotic + hetero, data = keeley)

# Use the `psem` function to create the SEM
model <- psem(abiotic_model, hetero_model, richness_model)

# Look at object
model

# Step 1: conduct tests of directed separation
# Establish the basis set & evaluate independence claims
# Missing path #1:
dsep1 <- lm(abiotic ~ hetero + distance, data = keeley)
# Missing path #2:
dsep2 <- lm(rich ~ distance + abiotic + hetero, data = keeley)

# Get P-values
P1 <- summary(dsep1)$coefficients[2, "Pr(>|t|)"]
P2 <- summary(dsep2)$coefficients[2, "Pr(>|t|)"]

# Construct C-statistic
C <- -2 * (log(P1) + log(P2))

C

# Compare to chi-squared distribution with 2*2 degrees of freedom
1 - pchisq(C, 4) # P < 0.05 == poor fit!

# Can use `dsep` function to perform the tests automagically
dSep(model)

# Can use `fisherC` function to evaluate claims
fisherC(model)

# The relationship between rich and distance is significant
# Re-introduce to the model
model2 <- update(model, rich ~ abiotic + hetero + distance)

model2

dSep(model2) # only 1 claim now
fisherC(model2) # P > 0.05 == model fits well!

# Get coefficients from good-fitting SEM
coefs(model2)

# Standardized estimates are in units of standard deviations of the mean
# Can be directly compared even though initial units are very different

# Plot SEM with standardized coefficients
plot(model2)

# Use `summary` function to get all information at once
summary(model2)

##### Chi-squared goodness of fit

# Evaluate fit using Chi-squared

# Sum Log-likelihoods from original model
LL_sem <- sum(sapply(model2[-4], logLik))

# Fit saturated model
saturated_model <- psem(
  lm(abiotic ~ distance + hetero, data = keeley),
  lm(hetero ~ distance, data = keeley),
  lm(rich ~ abiotic + hetero + distance, data = keeley)
)

# Sum Log-likelihoods from saturated model
LL_sat <- sum(sapply(saturated_model[-4], logLik))

Chi_sq <- -2 * (LL_sem - LL_sat)

Chi_sq

# Compare to chi-squared distribution with 1 df (one additional estimated 
# parameter in saturated model)
1 - pchisq(Chi_sq, 1) # P > 0.05 == good fit!

fisherC(model2) # slightly different than P-value from Fisher's C

# Re-fit in lavaan
library(lavaan)

form <- '
abiotic ~ distance 
hetero ~ distance
rich ~ abiotic + hetero + distance
'

sem(form, keeley) # same P-value!

# Can we test whether the model with the distance -> rich path 
# is statistically better?
AIC(model, model2)

anova(model, model2) # Chi-square difference test

#####-----Extensions to non-linear models---------------------------------------

# Re-fit Keeley example with GLM
model3 <- psem(
  lm(abiotic ~ distance, data = keeley),
  lm(hetero ~ distance, data = keeley),
  glm(rich ~ abiotic + hetero + distance, family = poisson(link = "log"), 
      data = keeley)
)

# Get summary
summary(model3)

# Compare with SEM of just LM 
anova(model2, model3) # GLM actually is less likely!

# Imagine that distance -> hetero relationship is truly nonlinear
# Re-fit SEM using generalized additive model (GAM)
library(mgcv)

model4 <- psem(
  lm(abiotic ~ distance, data = keeley),
  gam(hetero ~ s(distance), data = keeley),
  glm(rich ~ abiotic + hetero + distance, family = poisson(link = "log"), 
      data = keeley)
)

# Get summary 
summary(model4)

# Formally compare
anova(model2, model4) # Again, linear model is best
