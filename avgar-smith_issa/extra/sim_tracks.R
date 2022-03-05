################################################################X
#-------------ESA Statistical Methods Seminar Series------------X
#------Avgar & Smith -- Integrated Step Selection Analysis------X
#-------------------------07 March 2022-------------------------X
################################################################X

# Extra script: simulating movement data
# Brian J. Smith <brian.smith@usu.edu>

# Load packages ----
library(tidyverse)
library(raster)
library(amt)
library(lubridate)
library(circular)

# Habitat variables ----
hab <- stack("data/habitat.tif")
names(hab) <- c("forage", "temp", "predator", "cover")
plot(hab)

# ... movement-free habitat kernel ----
# We have to choose betas that will give us our movement-free habitat
# selection kernel.

# We interpret each beta as the log-RSS for a one-unit change in the covariate. 
# Note that the ratio of densities is the density of *steps* ending in that 
# habitat, GIVEN THE START LOCATION.

beta_forage = log(5)/500
beta_pred = log(0.25)/5
beta_temp2 = -1 * log(2)/36
beta_temp = beta_temp2 * -20

# Now let's add a bit of complexity. Let's say that our animal is nocturnal,
# and they spend the night foraging in grasslands and the day resting in
# forests.

# Let's stick with grassland as the reference category as it was in module 05.
# Specifically, we'll use grassland during daytime as the reference.
# We'll pick some betas such that:
#   Daytime:
#     Forest 10x grassland
#     Wetland 1/3 x grassland
#   Nighttime:
#     Forest 1/3x grassland 
#     Wetland 1/5x grassland 

# To be clear, we want the time of day at the end of the step to affect the
# habitat selected at the end of the step.

beta_forest_day = log(10)
beta_wetland_day = log(1/3)
beta_forest_night = log(1/3)
beta_wetland_night = log(1/5)

# ... selection-free movement kernel ----

# Now we have to specify the movement part of the model. 

# Let's assume our steps come from a gamma distribution and our turn
# angles come from a von Mises distribution.

# As we said above, our animal is nocturnal, so we want a step-length 
# distribution with short steps during the day and long steps at night.
# I.e., steps that start during day are short and steps that start during
# night are long.

# Let's choose these shape and scale parameters.
shp_day <- 2
scl_day <- 25

shp_night <- 5
scl_night <- 50

# Note that the mean of the gamma is given by shape * scale
# Day
shp_day * scl_day
# Night
shp_night * scl_night

# Plot the distributions
expand.grid(sl = seq(0.1, 1000, length.out = 100),
            time = c("day", "night")) %>% 
  mutate(
    shp = case_when(
      time == "day" ~ shp_day,
      time == "night" ~ shp_night),
    scl = case_when(
      time == "day" ~ scl_day,
      time == "night" ~ scl_night
    ),
    y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y, color = time)) +
  geom_line(size = 1) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()

# We also need to specify our turn angle distribution. Let's assume steps are
# more directed during at night and more uniform during the day.

# Von Mises parameters (assume mu = 0)
kappa_day <- 0.2
kappa_night <- 1.0

# Plot the distributions
expand.grid(ta = seq(-pi, pi, length.out = 100),
            time = c("day", "night")) %>% 
  mutate(
    kappa = case_when(
      time == "day" ~ kappa_day,
      time == "night" ~ kappa_night
    )) %>% 
# circular::dvonmises() is not vectorized
  rowwise() %>% 
  mutate(y = circular::dvonmises(ta, mu = 0, kappa = kappa)) %>% 
  ungroup() %>% 
  ggplot(aes(x = ta, y = y, color = time)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme_bw()

# ... times ----

# Now we need to define what time step our movement parameters correspond
# to. Let's say we have a GPS fix every hour and we want one month of data. 
# We'll straddle the spring equinox to have approximately the same number of 
# day and night hours to work with.

dat <- data.frame(id = "A01",
                  time = seq(ymd_hms("2021-03-05 0:00:00", tz = "US/Mountain"),
                             ymd_hms("2021-04-05 00:00:00", tz = "US/Mountain"),
                             by = "1 hour"),
                  x = NA,
                  y = NA)

# We'll start our animal just south of the middle of our map
dat$x <- mean(c(extent(hab)@xmin,extent(hab)@xmax))
dat$y <- mean(c(extent(hab)@ymin, extent(hab)@ymax))

# Check
plot(hab[[1]])
points(dat$x, dat$y, pch = 16, col = "red")

# We also want to assign the time of day to our locations. 
# `amt` can help with that.
?time_of_day

dat <- dat %>% 
  make_track(x, y, time, id = id, crs = 32612) %>% 
  time_of_day() %>% 
  # Back to a regular data.frame
  as.data.frame() %>% 
  # Rearrange/rename columns
  dplyr::select(id, x1 = x_, y1 = y_, t1 = t_, tod_start = tod_) %>% 
  # Convert to steps
  mutate(x2 = lead(x1),
         y2 = lead(y1),
         t2 = lead(t1),
         tod_end = lead(tod_start)) %>% 
  filter(!is.na(t2))

# Let's define our first step's endpoint as being 50m directly
# north to get us started moving. This is also the second step's start point.
dat$y2[1] <- dat$y1[2] <- dat$y1[1] + 50

# Lastly, let's add the absolute angle of the first step
dat$abs_angle <- NA
dat$abs_angle[1] <- 0 # directly north

# ... simulate movement ----
# Function to jitter data
jitter <- function(x, y, min = -25, max = 25) {
  res <- data.frame(x = x + runif(1, min, max),
                    y = y + runif(1, min, max))
  return(res)
}

# Function to crop raster to local coordinates
#   r = raster to crop
#   cent = vector of length 2 with xy-coordinates of centroid
#   d = distance (radius) of crop
crop_raster <- function(r, cent, d) {
  ext <- extent(cent[1] - d, cent[1] + d, cent[2] - d, cent[2] + d)
  res <- crop(r, ext)
  return(res)
}

# Maximum distance to calculate
qgamma(0.999, shape = shp_night, scale = scl_night)
max_d <- 800

# Loop over steps ----

# We already have the first step. We need to simulate the rest.
set.seed(20220307)

for (i in 2:nrow(dat)) {
  # Report status
  cat("\nStep", i, "of", nrow(dat))
  
  ## Is our step day or night?
  tod_st <- dat$tod_start[i]
  tod_end <- dat$tod_end[i]
  
  ## Calculate selection-free movement kernel
  # Start point
  start <- cbind(dat$x1[i], dat$y1[i])
  # Crop raster to 1000 m
  cropped <- crop_raster(hab, start, max_d)
  # Get coordinates of cropped raster
  coords <- coordinates(cropped)
  # Distances along x and y to every cell
  dx <- coords[, 1] - start[, 1]
  dy <- coords[, 2] - start[, 2]
  # Distance to every cell
  dists <- sqrt(dx^2 + dy^2)
  # Truncate at max distance
  trunc <- which(dists > max_d)
  dists[trunc] <- NA
  # Absolute angle to every cell
  abs <- (pi/2 - atan2(dy, dx)) %% (2*pi)
  # Relative angle difference
  rel_diff <- (abs - dat$abs_angle[i-1])
  # Relative angle
  rel_angle <- ifelse(rel_diff > pi, rel_diff - 2*pi, rel_diff)
  # Likelihood of step length
  sl_like <- dgamma(dists, 
                    shape = get(paste0("shp_", tod_st)),
                    scale = get(paste0("scl_", tod_st))) / (2 * pi * dists)
  # Necessary?
  sl_like <- sl_like/sum(sl_like, na.rm = TRUE)
  # Likelihood of turn angle (not vectorized -- need loop)
  ta_like <- rep(NA, length(sl_like))
  for (j in 1:length(ta_like)) {
    suppressWarnings({
      ta_like[j] <- dvonmises(rel_angle[j],
                              mu = 0,
                              kappa = get(paste0("kappa_", tod_st)))
      
    })
  }
  ta_like[trunc] <- NA
  # Necessary?
  ta_like <- ta_like/sum(ta_like, na.rm = TRUE)
  
  # Calculate kernel values
  move_kern_vals <- sl_like * ta_like
  # Normalize (sum to 1)
  move_kern_vals <- move_kern_vals/sum(move_kern_vals, na.rm = TRUE)
  
  # # If you want to plot this
  # move_kern <- cropped[[1]]
  # values(move_kern) <- move_kern_vals
  # plot(move_kern, main = "Movement Kernel")
  
  ## Movement-free habitat kernel
  hab_kern_vals <- as.data.frame(cropped, xy = TRUE) %>% 
    mutate(w = exp(
      beta_forage * forage +
        beta_temp * temp +
        beta_temp2 * temp^2 +
        beta_pred * predator +
        get(paste0("beta_forest_", tod_end)) * (cover == 2) +
        get(paste0("beta_wetland_", tod_end)) * (cover == 3)
    )) %>% 
    # Normalize
    mutate(w_prime = w/sum(w)) %>% 
    pull(w_prime)
  
  # # If you want to plot this
  # hab_kern <- cropped[[1]]
  # values(hab_kern) <- hab_kern_vals
  # plot(hab_kern, main = "Habitat Kernel")
  
  ## Combine
  step_kern_vals <- move_kern_vals * hab_kern_vals
  # Normalize
  step_kern_vals <- step_kern_vals/sum(step_kern_vals, na.rm = TRUE)
  
  # # If you want to visualize
  # step_kern <- cropped[[1]]
  # values(step_kern) <- step_kern_vals
  # plot(step_kern, main = "Habitat x Movement Kernel")
  
  # Randomly select cell to move into based on the probabilities
  cells <- 1:ncell(cropped)
  cells[trunc] <- NA
  next_cell <- sample(x = na.omit(cells),
                      size = 1,
                      prob = na.omit(step_kern_vals))
  
  # Get cell coordinates
  next_cell_coords <- xyFromCell(cropped, next_cell)
  
  # If you want to plot
  # points(next_cell_coords[,"x"], next_cell_coords[,"y"], pch = 16)
  
  # Jitter
  next_coords <- jitter(next_cell_coords[, 1], next_cell_coords[, 2])
  
  # Insert into data.frame
  if (i != nrow(dat)) {
    dat$x2[i] <- dat$x1[i+1] <- next_coords[, 1]
    dat$y2[i] <- dat$y1[i+1] <- next_coords[, 2]
  } else {
    dat$x2[i] <- next_coords[, 1]
    dat$y2[i] <- next_coords[, 2]
  }
  
  # Calculate absolute angle
  dx <- dat$x2[i] - dat$x1[i]
  dy <- dat$y2[i] - dat$y1[i]
  dat$abs_angle[i] = (pi/2 - atan2(dy, dx)) %% (2*pi)
  
  # If you want to check
  # dat[i, ]
  
}

# Save results
# dir.create("temp", showWarnings = FALSE)
# saveRDS(dat, "temp/sim.rds")

# Load results
# dat <- readRDS("temp/sim.rds")

# Now that we've simulated our data, let's see what our trajectory looks like.
traj <- dat %>% 
  filter(!is.na(t2)) %>% 
  select(id, x = x2, y = y2, time = t2) %>% 
  make_track(x, y, time, id = id, crs = 32612)

plot(hab[[1]])
points(traj)
lines(traj)

# Rename and save as CSV
traj2 <- traj %>%
  select(id, utm_e = x_, utm_n = y_, timestamp = t_) %>% 
  as.data.frame()

write.csv(traj2, "data/gps.csv", row.names = FALSE)
