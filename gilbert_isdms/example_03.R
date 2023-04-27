# Author: Neil Gilbert
# Date: 23 March 2023
# This script simulates abundance abundance of some focal critter across a landscape
# We then simulate distance sampling and count data from a portion of pixels in the landscape
# Then fit an integrated model combining the distance-sampling and count data
# You'll need the R packages listed below
# in addition, nimble requires a working C++ compiler - for example, Rtools on windows
# see http://r-nimble.org/html_manual/cha-installing-nimble.html

library(sf)
library(raster)
library(tidyverse)
library(fields)
library(patchwork)
library(nimble)
library(MCMCvis)

# number of pixels per edge of square landscape to simulate
npixels <- 50

# create a 50x50 (2500 cells) raster
# convert it to an sf polygons
# giving it a WGS84 CRS; this doesn't really matter
landscape <-
  raster::raster( nrows = npixels, 
                  ncols = npixels, 
                  xmn = -1, 
                  xmx =  1,
                  ymn = -1, 
                  ymx =  1 ) %>% 
  raster::rasterToPolygons(.) %>% 
  sf::st_as_sf(., crs = 4326, agr = "constant") %>% 
  dplyr::mutate(id = dplyr::row_number()) %>% 
  dplyr::select(id, geometry)

# calculate centroids of each pixel and
# add the x- and y- coordinates as columns
pixel_centroids <- landscape %>% 
  sf::st_centroid() %>% 
  dplyr::mutate(xcoord = sf::st_coordinates(.)[,1],
                ycoord = sf::st_coordinates(.)[,2])

#### Generate a continuous landscape pattern with splines ####
# The goal is to create a spatial covariate that we think
# will be associated with abundance of the simulated species
# think of this as % forest cover or terrain ruggedness

# create a systematic grid of knots across the landscape
knots <- base::expand.grid(
  x = base::seq( from = base::min( pixel_centroids$xcoord ) + 0.2,
                 to = base::max( pixel_centroids$xcoord ) - 0.2,
                 length.out = 4 ),
  y = base::seq( from = base::min( pixel_centroids$ycoord ) + 0.2,
                 to = base::max( pixel_centroids$ycoord ) - 0.2,
                 length.out = 4) )

# compute euclidean distances between pairs of points
omega_all <- ( fields::rdist( knots, knots ) ) ^ 3

svd_omega_all <- base::svd( omega_all )

sqrt_omega_all <- base::t( svd_omega_all$v %*% ( base::t( svd_omega_all$u ) * base::sqrt( svd_omega_all$d ) ) )

z_k <- ( fields::rdist( base::cbind(pixel_centroids$xcoord,
                                    pixel_centroids$ycoord),
                        knots) ) ^ 3

z_matrix <- base::t( base::solve( sqrt_omega_all, base::t( z_k ) ) )

z_matrix_scale <- base::apply(z_matrix, 2, base::scale)

# vector of spline coefficients
set.seed(11)
gamma <- stats::rnorm( base::nrow( knots ), 0, 0.5)

covariate <- z_matrix_scale %*% gamma

# visualize the covariate
landscape %>% 
  add_column(covariate = covariate) %>% 
  ggplot(aes(geometry = geometry, fill = covariate)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c()

#### simulate latent abundance ####
# Here, we will simulate the true, latent abundance of a species across the landscape

# abundance intercept
beta0 <- 1

# effect of the covariate
beta1 <- 0.5

# we're going to add in some noisiness to be more ecologically realistic
# this is the standard deviation of the normally distributed, zero-mean noise that we'll add
noisiness <- 0.2

# percentage of cells are sampled with distance sampling
prop_ds <- 0.03
# percentage of cells sampled with counts
prop_count <- 0.15

# add latent abundance to each cell as a function of the intercept, covariate and its coefficient, and noise
# also, add binary variables indicating whether the pixel was surveyed with distance sampling and counts 
true_abundance <- landscape %>% 
  tibble::add_column( covariate = covariate ) %>%
  dplyr::mutate( noise = stats::rnorm( n = base::nrow(.),
                                       mean = 0,
                                       sd = noisiness ), 
                 true_n = stats::rpois( n = base::nrow(.),
                                        lambda = base::exp( beta0 + beta1 * covariate + noise ) ) ) %>% 
  dplyr::mutate( ds_surveyed = stats::rbinom(n = base::nrow(.), size = 1, prob = prop_ds ), 
                 count_surveyed = stats::rbinom(n = base::nrow(.), size = 1, prob = prop_count ) ) 

# now, we need to create an expanded dataframe so that every individual animal has a row
# This is so we can assign distances to each animal
id_expanded <- c()
for( i in 1 : base::nrow( true_abundance ) ) {
  id_expanded <- c( id_expanded, base::rep( x = true_abundance[[i, "id"]], times = true_abundance[[i, "true_n"]] ) )
}

# maximum distance to which to count animals
# the value is arbitrary since this is a simulation
b <- 1
# how many distance bins to divide the surveyed area into
nbreaks <- 20
# width of each distance bin
width <- b / nbreaks

# intercept of function for scale parameter sigma,
# which controls how quickly detectabitily decays over distance
gamma0 <- -1
sigma <- base::exp( gamma0 )

expanded_df <- tibble::tibble( id = id_expanded ) %>% 
  # add distance of each individual to the transect
  # distances are uniformly distributed from on the transect (0) to the maximum distance to which animals are surveyed
  tibble::add_column( d = stats::runif(n = base::nrow(.), min = 0, max = b ) ) %>% 
  # calculate detection probability of individual based on its distance
  # we're using a half-normal detection function
  # also create dclass variable, which discretizes the distance measurements into bins
  dplyr::mutate( p = base::exp( -d * d / ( 2 * sigma * sigma ) ),
                 dclass = ( d %/% width ) + 1 ) %>%
  # create a column indicating whether the indivdiual was detected via distance sampling
  dplyr::mutate( ds_detected = stats::rbinom(n = nrow(.), size = 1, prob = p ) ) 

# breakpoints of the distance classes
dist_breaks <- base::seq( from = 0, to = b, by = width )

# the half-normal detection function
g <- function(x, sig) exp(-x * x / (2 * sig * sig))

# we'll calculate detection probability within each distance bin
p <- base::array( NA, dim = c( base::length( dist_breaks ) - 1 ) )
for(j in 1:base::length( p ) ) {
  p[j] <- ( stats::integrate( g,
                              dist_breaks[j],
                              dist_breaks[j + 1],
                              sig = sigma )$value / ( dist_breaks[j + 1] - dist_breaks[j] ) ) * ( width / b)
}

# the sum of bin-level probabilities can be considered overall detection probablity
# we'll use this to simulate the count data
p_counts <- sum( p )

survey_data <- expanded_df %>% 
  dplyr::group_by(id) %>% 
  # how many individuals were detected in the pixel with distance sampling?
  dplyr::summarise( n_detected_ds = sum( ds_detected ) ) %>% 
  dplyr::right_join( true_abundance ) %>% 
  # here, we simulate the pixel counts (without distance info) using the overall detection probability
  # we also convert to NA the counts for pixels that were not surveyed
  dplyr::mutate( count = stats::rbinom( n = base::nrow(.), size = true_n, prob = p_counts),
                 covariate = base::as.numeric( covariate ),
                 n_detected_ds = base::ifelse( ds_surveyed == 1, n_detected_ds, NA ), 
                 count = ifelse( count_surveyed == 1, count, NA ) ) %>% 
  dplyr::arrange(id) %>% 
  dplyr::select(
    id,
    covariate, 
    true_n, 
    ds_true = ds_surveyed,
    n_ds = n_detected_ds,
    count_true = count_surveyed,
    n_count = count,
    geometry = geometry) 

# let's visualize the spatial pattern of true abundance, as well as the two data sources

latent_n <- ggplot( survey_data, 
                    aes( geometry = geometry, fill = true_n ) ) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c() + 
  labs(title = "Latent abundance")+
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ds_plot <- ggplot(survey_data, 
                  aes( geometry = geometry, fill = n_ds ) ) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c(na.value = "white") + 
  labs(title = "Distance sampling: number detected")+
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

count_plot <- ggplot(survey_data, 
                     aes( geometry = geometry, fill = n_count ) ) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c(na.value = "white") + 
  labs(title = "Counts: number detected")+
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# use patchwork package to put the three plots together
latent_n | ds_plot | count_plot 

# package up the data and constants for the Nimble model
data <- list(
  # midpoint of the distance categories
  MIDPOINT = ( dist_breaks + width / 2 )[1:nbreaks],
  # width of each distance bin
  V = width, 
  # maximum distance to which animals are counted
  B = b,
  # distance classes of each individual detected with DS 
  DCLASS = dplyr::pull( dplyr::filter( expanded_df, ds_detected == 1), dclass),
  yDS = dplyr::pull( dplyr::filter( survey_data, !is.na(n_ds)), n_ds),
  yCount = dplyr::pull( dplyr::filter( survey_data, !is.na(n_count)), n_count),
  covariate = survey_data$covariate )

# constants for the model
# these are values used for indexing, looping, etc.
constants <- list(
  # number of distance bins
  NBINS = nbreaks,
  # the pixel index of cells surveyed with distance sampling
  ds_index  = dplyr::pull( dplyr::filter( survey_data, !is.na(n_ds)), id),
  # pixel index of cells surveyed with counts
  count_index = dplyr::pull( dplyr::filter(survey_data, !is.na(n_count)), id),
  # number of distance observations
  NDISTANCES = base::length(data$DCLASS),
  # number of pixels surveyed with distance sampling
  NDS = base::length(data$yDS),
  # number of pixels surveyed with counts
  NC = base::length(data$yCount)) 

# Now, here's the actual modeling part!
# The model is written in Nimble, which is a different language from R
# One thing to note is that this is a declarative language, meaning that
#--- order doesn't matter (unlike in an R script) 
# Nimble takes this code and compiles it to C++
code <- nimble::nimbleCode({
  
  ## PRIORS
  # using standard weakly informative priors
  beta0 ~ dnorm(0, sd = 2)  # abundance intercept
  beta1 ~ dnorm(0, sd = 2)  # effect of the covariate
  gamma0 ~ dnorm(0, sd = 2) # intercept for scale parameter in detection function
  
  ## DISTANCE SAMPLING SUBMODEL
  # first, the detection function
  sum_pie <- sum( pie[1:NBINS] ) # transect-level detection probability
  sigma <- exp( gamma0 ) # scale parameter in detection function - governs how detectability decays over distance
  
  # Loop through distance bins
  for( k in 1:NBINS ) { 
    log( g[k] ) <- -MIDPOINT[k] * MIDPOINT[k] / ( 2 * sigma * sigma ) # half-normal detection function
    pie[k] <- g[k] * ( V / B ) # multiply by the proportion of the whole surveyed zone that each bin represents
    pie_cell[k] <- pie[k] / sum_pie # compute cell probabilities
  }
  
  # Loop through the distance measurements
  for( i in 1:NDISTANCES ) {
    DCLASS[i] ~ dcat( pie_cell[1:NBINS] )
  }
  
  # Loop through pixels with distance-sampling data
  for( i in 1:NDS ){
    yDS[i] ~ dbin( sum_pie, N[i] )
    N[i] ~ dpois( lambda[i] )
    log( lambda[i] ) <- beta0 + beta1 * covariate[ ds_index[i] ] # model for abundance
  }
  
  ## COUNT SUBMODEL
  # loop through pixels with count data
  for( i in 1:NC ){
    yCount[i] ~ dbin( sum_pie, N_COUNT[i] )
    N_COUNT[i] ~ dpois( lambda_count[i] )
    log(lambda_count[i]) <- beta0 + beta1 * covariate[ count_index[i]] # note that beta0 and beta1 appear here again!
  }
})

# parameters to monitor
params <- c(
  "beta0",
  "beta1", 
  "sum_pie", 
  "gamma0")

# we have to provide initial values for the MCMC algorithm for 
# all the stochastic variables - i.e. those with "twiddles" in the code above that aren't data
inits <- function(){
  list(
    beta0 = rnorm(1, 0, 0.2), 
    beta1 = rnorm(1, 0, 0.2),
    gamma0 = runif(1, 0, 5), 
    N = data$yDS + 1,
    N_COUNT = data$yCount + 1)
}

# MCMC settings
nc <- 2
nburn <- 2000  
ni <- nburn + 3000
nt <- 3

# run the model!
out <- nimble::nimbleMCMC(
  code = code, 
  constants = constants, 
  data = data, 
  inits = inits(),
  monitors = params, 
  niter = ni, 
  nburnin = nburn, 
  nchains = nc, 
  thin = nt )

# inspect the output
# chains should be converged (Rhat < 1.1)
# 95% credible intervals should contain the true value
MCMCvis::MCMCsummary(out) %>% 
  tibble::as_tibble(rownames = "parameter") %>% 
  tibble::add_column( truth = c( beta0, beta1, gamma0, p_counts ) ) %>% 
  dplyr::select(parameter, truth, mean:n.eff)

# Let's map the model's predictions for abundance of the landscape, and compare total abundance (true) to estimate
# cross-join the covariate dataframe with posterior (cheating here and using only 1 chain)
all <- merge( dplyr::select( sf::st_drop_geometry( true_abundance ),
                             id, covariate),
              tibble::as_tibble( out[[1]][ ,1:2] ), # 1000 posterior samples of beta0 and beta1
              by = NULL)

# predicted lambda (i.e, expected abundance) for all pixels in landscape
predicted_lambda <- all %>% 
  dplyr::mutate( covariate = base::as.numeric( covariate ) ) %>% 
  dplyr::group_by( id ) %>% 
  dplyr::mutate( lambda = base::exp( beta0 + beta1 * covariate ) ) %>% 
  dplyr::summarise( lambda_mean = base::mean(lambda), 
                    lambda_l95 = stats::quantile(lambda, c(0.025)), 
                    lambda_u95 = stats::quantile(lambda, c(0.975)))

true_abundance %>% 
  dplyr::full_join( predicted_lambda ) %>% 
  ggplot(aes(geometry = geometry, fill  = lambda_mean)) + 
  geom_sf(color = NA) + 
  scale_fill_viridis_c()

# total abundance of landscape - hopefully truth should be contained in 95% CI
tibble::tibble(
  truth = base::sum( true_abundance$true_n ),
  totalN_l95 = base::sum( predicted_lambda$lambda_l95 ),
  totalN_mean = base::sum( predicted_lambda$lambda_mean ),
  totalN_u95 = base::sum( predicted_lambda$lambda_u95 ) )
