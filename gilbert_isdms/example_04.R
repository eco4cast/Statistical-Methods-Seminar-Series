# Integrated species distribution model (ISDM) tutorial

# In this tutorial, we will combine two datastreams within one model
# 1) Detection-nondetection data collected at fine spatial resolution (camera traps) 
# 2) Count data collected at coarse spatial resolution (harvest counts within counties)

# IMPORTANT - this script uses the NIMBLE R package
# This package is for Bayesian models; the models are compiled to C++ for speed
# Thus, you must have a working compiler installed on your machine in addition to the packages below
# Please  see https://r-nimble.org/html_manual/cha-installing-nimble.html for instructions

# Install these packages if you don't have them
# install.packages(c("raster", "tidyverse", "sf", "nimble"))

library(raster)
library(tidyverse)
library(sf)
library(nimble)

# set seed for entire session
addTaskCallback(function(...) {set.seed(1);TRUE})

# set resolution of the fine grid
# hereafter, "_f" suffix will refer to fine-resolution things 
nrow_f <- 36
ncol_f <- 36

# set resolution of the coarse grid
# hereafter, "_c" suffix will refer to coarse-resolution things 
nrow_c <- sqrt(nrow_f)
ncol_c <- sqrt(ncol_f)

# per-individual detection probability by cameras within fine-resolution cells
p_f <- 0.1

# per-individual harvest probability will vary randomly by coarse cell 
# each cell has probability of harvest somewhere in [0, 0.5]
# this mimics reality, since all individuals in a population are (hopefully) not harvested 
p_c <- runif(n = nrow_c*ncol_c, min = 0, max = 0.5)

# true probability of harvest is of course unknown
# instead, we rely on measures of harvest effort to inform probability of harvest
# assuming that increased effort will translate to higher probability of harvest
# we will say that we have some imperfect measure of harvest effort
# e.g., could be number of harvest authorizations or number of hunter hours in the field
noise <- p_c + rnorm(n = length(p_c), mean = 0, sd = 0.15)
effort <- (1 - 0.01)*((noise - min(noise))/(max(noise) - min(noise))) + 0.01
cor(p_c, effort) # correlation = 0.77 :)

# fine-resolution grid - blank raster, convert to poly for ease
fine <- raster(nrows = nrow_f, ncols = ncol_f, xmn = -1, xmx = 1, ymn = -1, ymx = 1) %>% 
  rasterToPolygons(fine) %>% 
  sf::st_as_sf(., crs = 4326, agr = "constant") %>% 
  mutate(id_f = row_number()) %>% 
  dplyr::select(id_f, geometry)

# coarse-resolution grid - blank raster, convert to poly for ease
coarse <- raster(nrows = nrow_c, ncols = ncol_c, xmn = -1, xmx = 1, ymn = -1, ymx = 1) %>% 
  rasterToPolygons(.) %>% 
  sf::st_as_sf(.,  crs = 4326, agr = "constant") %>% 
  mutate(id_c = row_number()) %>% # county id
  dplyr::select(id_c, geometry) %>% 
  add_column(p_c = p_c,
             effort = effort)

# what percentage of fine-resolution cells are NOT surveyed?
percent_na <- 0.85

#  intercept (b0) and slopes (b1 & b2) of linear function to generate abundance 
# IMPORTANT - these are the true values of parameters we are very interested in
# b0 is the intercept, or average expected abundance
# b1 and b2 are "species-environment relationships" aka coefficients for environmental predictors
# IRL, it is often of interest for ecologists managers to estimate these 
b0 <- 0.25
b1 <- -0.5
b2 <- 1.25

# identify rows that should NOT be surveyed
na_rows <- tibble(
  id_f    = 1:nrow(fine),
  cell_na = rbinom(nrow(fine), 1, percent_na))

# okay, big blob of code here...
# basically, what we're doing here is:
# 1) set up the true fine-res expected abundance values across the grid
# 2) generate 12 replicate surveys in a subset of the cells we're surveying with cameras
# 3) calculate coarse-cell abundance and generate harvest counts
# you'll get some warning messages, ignore em
(  final <-
    fine %>%  
    mutate(xcoord     = st_coordinates(st_centroid(.))[,1], # x-coordinate of cell
           ycoord     = st_coordinates(st_centroid(.))[,2], # y-coordinate of cell
           noise      = rnorm(n = nrow(.), mean = 0, 0.1),  # some cell specific noise
           log_lambda = b0 + b1*xcoord + b2*ycoord + noise, # abundance shall be greatest in upper left corner of grid
           lambda_f   = exp(log_lambda)) %>%                # expected abundance  
    st_drop_geometry() %>% 
      group_by(xcoord, ycoord) %>% 
    mutate(n_f        = rpois(n = 1, lambda = lambda_f),       # sample true abundance
           survey1    = rbinom(n = 1, size = n_f, prob = p_f), # generate twelve replicate surveys -
           survey2    = rbinom(n = 1, size = n_f, prob = p_f), # these surveys are binomial counts
           survey3    = rbinom(n = 1, size = n_f, prob = p_f), 
           survey4    = rbinom(n = 1, size = n_f, prob = p_f), 
           survey5    = rbinom(n = 1, size = n_f, prob = p_f), 
           survey6    = rbinom(n = 1, size = n_f, prob = p_f),
           survey7    = rbinom(n = 1, size = n_f, prob = p_f), 
           survey8    = rbinom(n = 1, size = n_f, prob = p_f), 
           survey9    = rbinom(n = 1, size = n_f, prob = p_f),
           survey10   = rbinom(n = 1, size = n_f, prob = p_f),
           survey11   = rbinom(n = 1, size = n_f, prob = p_f),
           survey12   = rbinom(n = 1, size = n_f, prob = p_f)) %>% 
    pivot_longer(survey1:survey12, names_to = "survey_name", values_to = "count") %>% 
    mutate(count = ifelse(count > 0, 1, 0)) %>% # "degrade" counts to detection/nondetection (1/0)
    full_join(na_rows) %>% 
    mutate(count = ifelse(cell_na == 1, NA, count)) %>%     # convert not-surveyed cells to NA
    ungroup() %>% 
    pivot_wider(names_from = "survey_name", values_from = "count") %>%
    ungroup(.) %>%   
    sf::st_as_sf(coords = c("xcoord", "ycoord"), crs = 4326, agr = "constant") %>%
    sf::st_join(coarse, join = st_within) %>% # we need to figure out which fine cells go in which coarse cells
    arrange(id_c) %>%
    group_by(id_c) %>%
    mutate(n_c = sum(n_f)) %>%  # coarse-cell abundance - sum of fine-cell abundance the cell contains
    mutate(harvest = rbinom(n = 1, size = n_c, p_c)) %>% # generate harvest count
    ungroup(.) %>% 
    # new ID - fine cells are no longer in "automatic" order b/c they are arranged by which coarse cell they fall in
    mutate(cell_id = row_number()) %>% 
    st_drop_geometry() %>% 
    full_join(fine) %>%
    mutate(xcoord     = st_coordinates(st_centroid(geometry))[,1], # x-coordinate of cell
           ycoord     = st_coordinates(st_centroid(geometry))[,2]) %>% 
    
    dplyr::select(cell_id, id_f, id_c,
                  xcoord, ycoord,
                  lambda_f, n_f,
                  survey1:survey12, n_c, p_c, effort, harvest, geometry) )

glimpse( final )

# Okay, let's visualize the simulated data

# first, this is fine-resolution EXPECTED abundance
# note, expected abundance is continuous (i.e. can be decimal number)
# should increase fairly smoothly to the upper left of the grid
final %>% 
  ggplot(aes(geometry = geometry, fill = lambda_f)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c("Expected abundance") + 
  theme_void() 

# now, TRUE abundance within fine cells
# looks choppier because of natural Poisson randomness
# also, note that true abundance is an integer (it's a count)
final %>% 
  ggplot(aes(geometry = geometry, fill = n_f)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c("Abundance") +
  theme_void() 

# now, detection-nondetection surveys within fine cells
# remember, only 15% of the cells are surveyed
final %>% 
  dplyr::select(id_f, geometry, survey1:survey12) %>% 
  pivot_longer(survey1:survey12, names_to = "survey", values_to = "detect") %>%
  mutate(survey = factor(survey, levels = c("survey1", "survey2", "survey3",
                                            "survey4", "survey5", "survey6", 
                                            "survey7", "survey8", "survey9",
                                            "survey10", "survey11", "survey12"))) %>% 
  ggplot(aes(geometry = geometry, fill = factor(detect))) + 
  geom_sf(color = NA) + 
  facet_wrap(~survey) +
  scale_fill_viridis_d("Detected?", 
                       labels = c("No", 
                                  "Yes", 
                                  "Not surveyed"),
                       option = "B", 
                       begin = 0.15, 
                       end = 0.85,
                       na.value = "gray90") + 
  theme_void() 

# coarse-resolution true abundance
final %>% 
  ggplot(aes(geometry = geometry, fill = n_c)) + 
  geom_sf(color = NA) +
  scale_fill_viridis_c("Abundance") +
  theme_void()

# harvest-probability heterogeneity
# this could be due to variation in hunting behavior by humans - areas w/
# higher p(harvest) might, for example, have more public land, providing more access to animals
final %>% 
  ggplot(aes(geometry = geometry, fill = p_c)) +
  geom_sf(color = NA) + 
  scale_fill_viridis_c("P(Harvest)",
                       option = "B", 
                       begin = 0.2, 
                       end = 0.8) +
  theme_void()

# and finally, our coarse-resolution harvest counts
# note how the harvest counts are an imperfect representation of the true abundance pattern
final %>% 
  ggplot(aes(geometry = geometry, fill = harvest)) + 
  geom_sf(color = NA) +
  scale_fill_viridis_c("Harvest") +
  theme_void()

# moving on...
# prepare the coarse-cell data
# id_c is the identity of the coarse cell
# effort is our imperfect representation of harvest effort 
# harvest is the harvest count
# low and high are indices to tell the model which coarse cell each fine cell falls inside
( coarse_indices <- final %>%
    dplyr::select(id_f, id_c, effort, harvest) %>%
    st_drop_geometry(.) %>%
    ungroup(.) %>%
    mutate(row = row_number()) %>%
    group_by(id_c) %>%
    mutate(low = min(row), 
           high = max(row)) %>%
    dplyr::select(-id_f, -row) %>%
    distinct(.) )

# prepare the fine-resolution data
# basically, we just want to grab the survey data for the cells that are indeed surveyed
( surveyed_cells <- final %>%
    filter(!is.na(survey1)) %>%
    ungroup(.) %>%
    dplyr::select(id_f, survey1:survey12, cell_id) %>%
    sf::st_drop_geometry(.) ) 

# we have to package up the data for the model
data <- list(
  xcoord  = as.vector(final$xcoord), # predictor
  ycoord  = as.vector(final$ycoord), # predictor
  y       = unname(as.matrix(dplyr::select(surveyed_cells, starts_with("survey")))), # detection-nondetection - response
  harvest = coarse_indices$harvest, # harvest count - response data
  effort  = coarse_indices$effort) # effort data

# similarly, we have to package up the "constants" for the model
# these are pretty much all used for indexing purposes within loops and such
constants <- list(
  ncell_f        = length(data$xcoord),    # number of fine-resolution cells
  ncell_surveyed = nrow(data$y),           # number of fine-resolution cells that are surveyed
  nsurveys       = ncol(data$y),           # number of replicate surveys
  cell           = surveyed_cells$cell_id, # the cell id's of the cells that are surveyed
  ncell_c        = nrow(coarse_indices),   # number of coarse cells
  low            = coarse_indices$low,     # lower index to say which fine cells are in a given coarse cell
  high           = coarse_indices$high)    # upper index to say which fine cells are in a given coarse cell

# now, the model 
# important to thing to note
# with the model, we are creating an object that is written in a different language (i.e., not R)
# it is nimble, which is a dialect of BUGS (Bayesian inference Using Gibbs Sampling)
code <- nimble::nimbleCode({
  
  # Priors
  # Basically, we're giving the model an idea of what the distribution of a parameter will look like
  # e.g., will it be in [0,1], unconstrained, etc
  b0     ~ dnorm(0, sd = 2) # intercept
  b1     ~ dnorm(0, sd = 2) # coefficient for x-coord
  b2     ~ dnorm(0, sd = 2) # coefficient for y-coord
  p      ~ dbeta(1, 1)      # detection probability - constraining this to [0, 1]
  gamma0 ~ dnorm(0, sd = 2) # intercept for abundance-scaling function
  gamma1 ~ dnorm(0, sd = 2) # slope for abundance-scaling function
  
  # Likelihood
  # SDM - loop through all fine-resolution cells
  for(i in 1:ncell_f){
    log(lambda[i]) <- b0 + b1*xcoord[i] + b2*ycoord[i]
    n[i] ~ dpois(lambda[i]) # latent true abundance
  }
  
  # camera submodel - detection prob is held constant
  # importantly, this model assumes that detection is conditional on abundance
  # so if a camera detects an animal every survey, that location probably has higher abundance
  for(j in 1:ncell_surveyed){
    for(k in 1:nsurveys){
      muy[j, k] <- 1 - pow(1 - p, n[cell[j]])
      y[j, k] ~ dbern(muy[j, k])
    }}
  
  # harvest coarse-resolution submodel
  for(q in 1:ncell_c){
    # function to scale fine-scale to county-scale expected abundance
    log(lambda_county[q]) <- gamma0 + gamma1*log(sum(lambda[low[q]:high[q]])) 
    harvest[q] ~ dpois(effort[q]*lambda_county[q])
  }
})

# MCMC settings - burnin and total number of iterations
# basically, how long do we want to run the model for 
nb <- 24000
ni <- nb + 1000

# parameters that we want to keep track of
keepers <- c("lambda", "b0", "b1", "b2")

# we have to give the parameters within the model initial values 
# TLDR; the algorithm needs somewhere to start
inits <- function() {
  base::list(n      = rep(1, constants$ncell_f),
             b0     = runif(1, -1, 1),
             b1     = runif(1, -1, 1),
             b2     = runif(1, -1, 1),
             gamma0 = runif(1, -1, 1), 
             gamma1 = runif(1, -1, 1),
             p      = runif(1, 0, 1))}

# run the model! 
# this will take a little bit
# it will take a minute or two for the model to compile (i.e. translate to C++)
# then, you'll see a progress bar pop up once that happens
# all told, will take a few minutes (depends on if you change the iterations), so grab a snack or something
out <- nimble::nimbleMCMC(code = code,
                          data = data, 
                          constants = constants, 
                          monitors = keepers,
                          inits = inits(), 
                          nburnin = nb, 
                          ni = ni)

# grab the results for the estimated fine-cell abundances
# important note - we get a posterior DISTRIBUTION for each lambda
# so, for simplicity, we'll grab those and summarise the median value
lambda_result <- as_tibble(out) %>%
  dplyr::select(starts_with("lambda")) %>%
  pivot_longer(1:constants$ncell_f, names_to = "cell", values_to = "lambda_est") %>%
  mutate(cell_id = as.integer(extract_numeric(cell))) %>%
  group_by(cell_id) %>%
  summarise(lambda_est = median(lambda_est)) %>%
  ungroup(.) %>% 
  left_join(dplyr::select(final, cell_id, lambda_actual = lambda_f, geometry) %>% ungroup(.))

# let's compare the model's estimates for the expected abundance against the true values
ggplot(lambda_result, aes(x = lambda_actual, y = lambda_est)) + 
  geom_point(alpha = 0.25) + 
  geom_abline(slope = 1, intercept = 0, color = "red")

# map it out - does the spatial pattern of expected abundance predicted by the model match reality?
lambda_result %>% 
  pivot_longer(lambda_est:lambda_actual, names_to = "which", values_to = "lambda") %>% 
  ggplot(aes(geometry = geometry, fill = lambda)) + 
  geom_sf(color = NA) + 
  facet_wrap(~which) + 
  scale_fill_viridis_c()

# same thing for the parameters b0, b1, b2
# recall, these would represent species-environment relationships that are typically of interest
# we'll summarise the mean and the 95% credible interval
betas <- as_tibble(out) %>% 
  dplyr::select(b0, b1, b2) %>% 
  pivot_longer(b0:b2, names_to = "parameter", values_to = "value") %>% 
  group_by(parameter) %>% 
  summarise(mean = mean(value), 
            lower = quantile(value, 0.025), 
            upper = quantile(value, 0.975))

# how'd the model do estimating those parameters?
# red vertical lines mark the true values of each parameter
ggplot(betas, aes(x = mean, y = parameter)) + 
  geom_errorbar(aes(xmin = lower, xmax = upper), size = 2, width = 0) + 
  geom_vline(xintercept = b0, color = "red") + 
  geom_vline(xintercept = b1, color = "red") + 
  geom_vline(xintercept = b2, color = "red")
