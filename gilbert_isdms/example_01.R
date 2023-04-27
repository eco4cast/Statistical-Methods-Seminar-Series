# Date: 20 April 2023
# Author: Neil Gilbert
# iSDM combining replicated counts and single-visit presence-absence data

library(spatstat)
library(tidyverse)
library(sf)
library(nimble)
library(MetBrewer)
library(MCMCvis)

# simulation an inhomogenous point pattern - increasing intenstiy of points with x coordinate
d <- spatstat.random::rpoispp( function(x, y) { 25 + 15 * x },
                               win=owin( c( -2, 2 ), c( -2, 2) ), nsim = 1 )

# put the points into an sf object
pts <- tibble(
  x = d$x,
  y = d$y) %>% 
  st_as_sf(., coords = c("x", "y"), crs = 4326)

npixels <- 25
# proportion of cells with count data
prop_counts <- 0.05
# proportion of cells with presence/absence data
prop_det <- 0.35
# detection probability
detection_prob <- 0.5

# create a 25x25 (625 cells) raster
# convert it to an sf polygons
# giving it a WGS84 CRS; this doesn't really matter
landscape <- raster::raster( nrows = npixels, 
                             ncols = npixels, 
                             xmn = -2, 
                             xmx =  2,
                             ymn = -2, 
                             ymx =  2 ) %>% 
  raster::rasterToPolygons(.) %>% 
  sf::st_as_sf(., crs = 4326, agr = "constant") %>% 
  dplyr::mutate(id = dplyr::row_number()) %>% 
  dplyr::select(id, geometry) %>% 
  # some acrobatics to add, for each cell, the number of points falling inside of it
  dplyr::full_join(
    sf::st_join( pts, ., join = st_within) %>% 
      sf::st_drop_geometry() %>% 
      dplyr::count(id)
  ) %>% 
  dplyr::mutate( n = ifelse(is.na(n), 0, n)) %>% 
  # grabbing the xcoordinate of the pixel centroid (to be used as a predictor)
  dplyr::mutate( xcoord = st_coordinates( st_centroid(.))[,1]) %>% 
  dplyr::rowwise() %>% 
  # is the cell surveyed with repeated counts?
  dplyr::mutate(surveyed_counts = rbinom( 1, 1, prop_counts)) %>% 
  # is the cell surveyed with single-visit detection/nondetection data?
  dplyr::mutate(surveyed_det = ifelse(surveyed_counts == 0, rbinom(1, 1, prop_det), 0)) 

ggplot() +
  geom_sf(data = landscape, aes(geometry = geometry)) +
  geom_sf(data = pts, aes(geometry = geometry))


count <- landscape %>% 
  # simulate count data
  dplyr::mutate( count = rbinom(1, size = n, prob = detection_prob)) %>% 
  # knock out data for unsurveyed cells
  dplyr::mutate( count = ifelse(surveyed_counts == 0, NA, count)) %>% 
  dplyr::filter( !is.na(count)) %>% 
  # generate replicate counts - assuming detection probability is 0.5
  dplyr::mutate( `Visit 1` = rbinom( 1, n, detection_prob),
                 `Visit 2` = rbinom( 1, n, detection_prob), 
                 `Visit 3` = rbinom( 1, n, detection_prob), 
                 `Visit 4` = rbinom( 1, n, detection_prob)) %>% 
  dplyr::select( id, n, `Visit 1`:`Visit 4` ) %>% 
  tidyr::pivot_longer(`Visit 1`:`Visit 4`, names_to = "rep", values_to = "count") 

# visualize the count data
( count_plot <- ggplot( ) +
    geom_sf(data = count, aes( geometry = geometry, fill = factor(count)), color = NA) +
    geom_sf(data = pts, aes(geometry = geometry), alpha = 0.2) +
    facet_wrap(~rep) + 
    theme_void() +
    scale_fill_manual("Count",
                      values = MetBrewer::MetPalettes$Benedictus[[1]][c(6,5,4,3,2,1)]) )

# simulate the presence/absence data
detection_nondetection <- landscape %>% 
  rowwise() %>% 
  mutate( z = ifelse( n > 0, 1, 0)) %>% 
  mutate( det = rbinom(1, z, prob = 0.5)) %>%  
  mutate( det = ifelse(surveyed_det == 0, NA, det)) %>% 
  filter( !is.na(det)) %>% 
  dplyr::select(id, n, z, det) 

( detection_nondetection_plot <- ggplot() +
    geom_sf( data = detection_nondetection, aes(geometry = geometry, fill = factor(det)), color = NA) +
    geom_sf(data = pts, aes(geometry = geometry)) +
    scale_fill_manual( "Detected", 
                       values = MetBrewer::MetPalettes$Hiroshige[[1]][c(3, 7)], 
                       labels = c("No", "Yes")) +
    theme_void() )

# format the count data for input in the model
count_wide <- count %>% 
  st_drop_geometry() %>% 
  pivot_wider(names_from = rep, values_from = count)

# package up the data to go into the model
data <- list(
  detection = detection_nondetection$det, # presence absence data
  count = unname ( as.matrix ( select(count_wide, 3:6 ) ) ), # replicated count data
  xcoord = unname(landscape$xcoord)) # xcoord of pixel centroid - to be used as a predictor

# package up model "constants" - mostly used for indexing
constants <- list(
  det_id = detection_nondetection$id, # id of cells with presence/absence data
  count_id = count_wide$id, # id of cells with count data
  npixels = length(data$xcoord), # total number of cells in landscape
  npixels_with_counts = nrow(data$count), # number of cells with count data
  nreps = ncol(data$count), # number of repeat visits for counts
  npixels_with_detection = length(data$detection)) # number of cells w/ presence/absence data

# here is the nimble model code
code <- nimble::nimbleCode({
  
  # Priors
  #-------------------------------
  beta0 ~ dnorm(0, sd = 5)
  beta1 ~ dnorm(0, sd = 2)
  p ~ dbeta(1, 1)
  
  # Model for ecological state (SDM)
  #--------------------------------------
  for( i in 1:npixels ) {
    log( lambda[i] ) <- beta0 + beta1 * xcoord[i]
    N[i] ~ dpois( lambda[i] )
  }
  
  # Observation submodel for the replicated counts
  #-----------------------------------------------
  for( i in 1:npixels_with_counts ) {
    for( j in 1:nreps ) {
      count[i, j] ~ dbin(p, N[count_id[i]])
    }
  }
  
  # Observation submodel for the presence-absence data
  #-----------------------------------------------------
  for( i in 1:npixels_with_detection ) {
    detection[i] ~ dbern( p * ( 1 - exp( - lambda[ det_id[i] ] ) ) )
  }
  
  # get total abundance within landscape as a derived variable
  totalN <- sum( N[1:npixels] )
  
})

# we have to give the MCMC algorithm starting values for parameters
# for the latent N, we have to max sure that the initial values don't
# conflict with the observed data - e.g., that the initial value is less
# than the observed count - so for cells with counts, put the initial value
# as the max + 1 of the observed counts
Nst <- rpois(constants$npixels, 1)
Nst[constants$count_id] <- apply(data$count, 1, max) + 1

inits <- function(){
  list(
    beta0 = rnorm(1, 0, 1), 
    beta1 = rnorm(1, 0, 1), 
    p = runif(1, 0.1, 0.9), 
    N = Nst
  )
}

# parameters we want to track
params <- c("beta0", "beta1", "p", "totalN")

# for the sake of speed here, we are going to run 1 very short MCMC chain
# in practice, you will want to run multiple chains, probably for much longer
out <- nimble::nimbleMCMC(code = code,
                          data = data, 
                          constants = constants, 
                          monitors = params,
                          inits = inits(), 
                          nburnin = 2000, 
                          ni = 5000,
                          nchains = 2)

MCMCvis::MCMCsummary(out, params = c("totalN", "p", "beta0", "beta1"))
