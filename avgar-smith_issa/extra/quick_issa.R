################################################################X
#-------------ESA Statistical Methods Seminar Series------------X
#------Avgar & Smith -- Integrated Step Selection Analysis------X
#-------------------------07 March 2022-------------------------X
################################################################X

# Main script: quick iSSA using piped workflow
# Brian J. Smith <brian.smith@usu.edu>

# Load packages ----
library(tidyverse)
library(raster)
library(amt)
library(lubridate)

# Load data ----
# Load GPS data
gps <- read_csv("data/gps.csv")
# Set timezone to correct value (US/Mountain)
tz(gps$timestamp) <-  "US/Mountain"

# Load habitat data
hab <- stack("data/habitat.tif")
# Give raster layers names
names(hab) <- c("forage", "temp", "predator", "cover")

# Fit iSSF ----
set.seed(123 * 456)
issf <- gps %>% 
  make_track(utm_e, utm_n, timestamp, id = id, crs = 32612) %>% 
  track_resample(rate = hours(1), tolerance = minutes(5)) %>% 
  steps_by_burst() %>% 
  random_steps(n_control = 20) %>% 
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_)) %>% 
  extract_covariates(hab) %>% 
  mutate(cover = factor(cover,
                     levels = 1:3,
                     labels = c("grassland", "forest", "wetland"))) %>% 
  time_of_day(where = "both") %>% 
  # Make your own dummy variables
  # Grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night")) %>%
  # Fit the model
  fit_issf(case_ ~ 
             # Habitat
             forage + temp + I(temp^2) + predator + 
             # All the cover terms
             forest_day + wetland_day + 
             forest_night + wetland_night +
             # Movement
             tod_start_ : log_sl_ + tod_start_ : sl_ + 
             cos_ta_ + tod_start_ : cos_ta_ +
             # Strata (steps)
             strata(step_id_), model = TRUE)

# Model summary
summary(issf)
