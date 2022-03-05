################################################################X
#-------------ESA Statistical Methods Seminar Series------------X
#------Avgar & Smith -- Integrated Step Selection Analysis------X
#-------------------------07 March 2022-------------------------X
################################################################X

# Main script: iSSA demo using 'amt'
# Brian J. Smith <brian.smith@usu.edu>

# In this demo, we will read in some simulated data and format it using
# 'amt'. Then we will address the steps to implementing iSSA:
#   1. Sample
#   2. Attribute
#   3. Analyze
#   4. Infer
#   5. Predict (not covered here)

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
# Quick look
print(gps, n = 5)

# Load habitat data
hab <- stack("data/habitat.tif")
# Give raster layers names
names(hab) <- c("forage", "temp", "predator", "cover")
# Quick look
plot(hab)

# Formatting data with 'amt' ----
# The basic building block in the package 'amt' is a 'track_*' object.
# For movement data, this is a 'track_xyt' object, but locations without
# timestamps are also possible ('track_xy').

# We can format our data as a 'track_*' object using 'make_track()'.
trk <- make_track(gps, utm_e, utm_n, timestamp, id = id, crs = 32612)

# We can see the class of the resulting object
class(trk)

# This is still a tibble and a data.frame, so we can manipulate it just
# like we would manipulate any other tibble or data.frame in R. But we
# get some added benefits from it being a 'track_xyt' object -- e.g.,
# we have default plotting methods:
plot(trk)
lines(trk)

# Or, on top of our raster:
plot(hab[[1]])
points(trk)
lines(trk)

# Turning tracks into steps ----
# While the most basic building block in 'amt' is a track, iSSA requires
# that we convert from a point representation to a step representation.

# We can do that with the function 'steps()'.
stp <- steps(trk)
# View(stp)

# Note that steps have their own S3 class:
class(stp)

# Data management ----
# Data cleaning is an important first step whenever you're working with
# GPS data. Unfortunately, we do not have time in this webinar to dive
# into the details of data cleaning. We have data cleaning functions in 
# 'amt', and we will soon have a new vignette with a demonstrated cleaning
# workflow.

# Once your data are cleaned, you will often have gaps in an otherwise
# regular trajectory. iSSA requires a constant step duration, so we cannot
# create steps across gaps in our data. Let's demonstrate how to handle
# that in 'amt'.

# Randomly remove ~5% of our locations:
set.seed(12345)
trk2 <- trk %>% 
  mutate(rm = as.logical(rbinom(n = nrow(trk), size = 1, prob = 0.05))) %>% 
  filter(!rm)

# Check
nrow(trk2)/nrow(trk)
print(trk2, n = 10)

# Notice we removed our 10th location -- we have a 2h gap from 
# 2021-03-05 09:00 to 2021-03-05 11:00. 

# What happens if we make steps, now?
stp2 <- steps(trk2)
# View(stp2)

# We now have some steps where dt = 2 hours. Our trajectory is irregular, 
# and that will not work for iSSA. We can filter out those 2 hour steps
# fairly easily in this simplified example, but a more powerful way to
# deal with this is to use the function 'track_resample()'.

# 'track_resample()' takes the track, the desired duration, and a tolerance
# around that duration as arguments, then it divides a trajectory into
# "bursts" with a constant step duration. Note, this function doesn't do
# any interpolation.
trk3 <- track_resample(trk2, rate = hours(1), tolerance = minutes(5))
print(trk3, n = 12)

# We have a new column, 'burst_', that identifies our bursts. We can see the
# switch from burst 1 to burst 2 after our gap at 10:00. If we pass that to
# 'steps()', we now get a warning:
steps(trk3)

# The warning tells us to use 'steps_by_burst()' to make use of the new 
# column:
stp3 <- steps_by_burst(trk3)
# View(stp3)

# We no longer have any steps with 2-h durations, and we are ready to
# proceed with our iSSA.

# Step 1. Sample ----
# The next thing we need to do is generate random available steps from
# a tentative parametric distribution. We will use the gamma distribution
# to generate step lengths and a von Mises distribution to generate turn
# angles. After fitting our iSSA, we will update these tentative 
# distributions to the true, selection-free movement distributions.

# A good way to decide on tentative parameters is to fit them to observed
# steps. We can do that in 'amt' like this:

# gamma distribution (fit to step lengths)
fit_distr(stp3$sl_, "gamma")

# von Mises distribution (fit to turn angles)
fit_distr(stp3$ta_, "vonmises")

# We can use the function 'random_steps()' to generate our random steps.
# When we pass a 'steps_xyt' object, the default is for it to fit the
# gamma and von Mises distributions, just as we did above.
set.seed(20220307)
obs_avail <- random_steps(stp3, n_control = 20)

# Random steps have their own S3 class, but they are also still of
# class 'steps_xyt'.
class(obs_avail)

# Take a look:
# View(obs_avail)

# Each observed step gets an identifying number (step_id_), and all of 
# the paired available steps get that same ID. These IDs form the strata
# in our conditional logistic regression. The variable 'case_' is TRUE
# for the observed step and 'FALSE' for the available steps. This will
# be the response variable in our conditional logistic regression.

# Note that the tentative step length and turn angle distributions are
# attached to the object as attributes.
attributes(obs_avail)

# You can also access them with these convenience functions:
sl_distr(obs_avail)
ta_distr(obs_avail)

# Since we used the gamma distribution as our tentative step-length
# distribution, we need to include step length and log(step length)
# in our iSSF. We already have a column, 'sl_', but we need to add 
# 'log_sl_'. Likewise, we need the cosine of the turn angle to update
# the von Mises distribution, so we need to add cos_ta_:
obs_avail <- obs_avail %>% 
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_))

# Have a look:
print(obs_avail, n = 3, width = 100)

# Step 2. Attribute ----

# ... attach environmental covariates ----
# Now that we have observed and available steps, we need to attach our
# environmental covariates to each one. We can use the functions
# 'extract_covariates()' to attach the raster values to our steps.
covs <- extract_covariates(obs_avail, hab)

# Note that 'cover' is a factor, and R doesn't know that right now.
# Let's format it properly.
covs$cover <- factor(covs$cover,
                     levels = 1:3,
                     labels = c("grassland", "forest", "wetland"))

# Have a look:
print(covs, n = 3, width = 100)

# ... attach time of day ----
# One of the great strengths of iSSA is that it can handle temporal 
# variation. Perhaps our organism selects habitat differently between
# day and night, or perhaps it moves with different speeds during day
# or night. We can account for this with interactions in our iSSF,
# and 'amt' has a function that extracts the time of day,
# given the coordinates and the date. I.e., it accounts for different
# sunset times at different latitudes on different days of the year.
covs2 <- time_of_day(covs, where = "both")

# The argument 'where = "both"' tells the function to extract the time
# of day for *both* the start of the step and the end of the step.
print(covs2, n = 3, width = 100)

# Now we're ready to fit a model!

# Step 3: Analyze ----
# Let's begin with a simple iSSF. We'll model our movement-free habitat
# selection kernel as a function of forage and predation risk, and we'll 
# include all the movement parameters to update the selection-free movement 
# kernel.

# Note that the 'amt' function 'fit_issf()' is just a wrapper for
# survival::clogit(). For much more information on the implementation,
# see the help file.
?survival::clogit

m1 <- fit_issf(covs2, 
               # Response
               case_ ~ 
                 # Habitat
                 forage + predator +
                 # Movement
                 sl_ + log_sl_ + cos_ta_ + 
                 # Stratum
                 strata(step_id_),
               # Need this later for model predictions
               model = TRUE)

# Let's have a look at the structure of our object.
str(m1, 1)

# Notice that it is a list with 4 elements at the top level.
#   - $model: the actual fitted model
#   - $sl_: the tentative step-length distribution
#   - $ta_: the tentative turn-angle distribution
#   - $more: (currently empty) a placeholder for additional information

# Step 4: Infer ----
# Now we're ready to draw some inference from our simple model.

# Take a look at the model summary:
summary(m1)

# Main take-aways from this summary are:
#   (1) our animal selects for forage (+ coef)
#   (2) our animal avoids predation risk (- coef)
#   (3) our tentative distributions are very close to the true selection-free
#       movement kernel (none of the movement params are significant.)

# ... selection-free movement kernel ----
# With this simple model structure, we have an 'amt' function that will
# automatically update the movement distributions.

# Tentative step-length distribution
(tent_sl <- sl_distr(m1))
# Updated selection-free step-length distribution
(upd_sl <- update_sl_distr(m1))

# Compile the parameters into a data.frame for plotting with ggplot
tent_df <- data.frame(dist = "tent",
                      shp = tent_sl$params$shape,
                      scl = tent_sl$params$scale)
upd_df <- data.frame(dist = "upd",
                     shp = upd_sl$params$shape,
                     scl = upd_sl$params$scale)

(sl_df <- rbind(tent_df, upd_df))

# Plot
expand.grid(sl = seq(1, 1000, length.out = 100),
            dist = c("tent", "upd")) %>% 
  left_join(sl_df) %>% 
  mutate(y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y, color = dist)) +
  geom_line() +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()

# We can see there is barely a difference between our tentative and updated
# step-length distributions. This is no surprise, given that the betas for
# sl_ and log_sl_ were not significant.

# Tentative turn-angle distribution
(tent_ta <- ta_distr(m1))
# Updated selection-free turn-angle distribution
(upd_ta <- update_ta_distr(m1))

# Compile the parameters into a data.frame for plotting with ggplot
tent_df_ta <- data.frame(dist = "tent",
                         k = tent_ta$params$kappa)
upd_df_ta <- data.frame(dist = "upd",
                        k = upd_ta$params$kappa)

(ta_df <- rbind(tent_df_ta, upd_df_ta))

# Plot
expand.grid(ta = seq(-pi, pi, length.out = 100),
            dist = c("tent", "upd")) %>% 
  left_join(ta_df) %>% 
  # circular::dvonmises is not vectorized
  rowwise() %>% 
  mutate(y = circular::dvonmises(ta, mu = 0, kappa = k)) %>% 
  ggplot(aes(x = ta, y = y, color = dist)) +
  geom_line() +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme_bw()

# Again, very little difference.

# ... movement-free selection kernel ----
# We can use (log-)RSS to visualize the habitat selection part of the model.
# We have a function in 'amt' to calculate log-RSS(x1, x2) for any x1 and
# x2 of interest.

# Both x1 and x2 are defined using data.frames and require all of the
# covariates in the fitted model. x1 can be any number of rows you wish,
# but to avoid unintended issues with R's recycling rules, we limit
# x2 to be exactly 1 row.

## Scenario 1:
# How much more likely is our animal to step into a habitat with 
# forage = 600 g/m^2 than a habitat with forage = 400 g/m^2?
x1 <- data.frame(forage = 600, predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

x2 <- data.frame(forage = 400, predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

logRSS <- log_rss(m1, x1 = x1, x2 = x2)

# RSS
exp(logRSS$df$log_rss)

# Our animal is ~ 1.5x more likely to step into a habitat with forage = 600
# than forage = 400.

# Note the structure of the object we created.
str(logRSS, 1)

# It is a list with 4 elements.
#  - 'df' -- a data.frame with the results, useful for plotting
#  - 'x1' -- the original x1 passed to the function
#  - 'x2' -- the original x2 passed to the function
#  - 'formula' -- the formula of the model

# The 'df' element is the one you'll typically be interested in.

## Scenario 2:
# How much more likely is our animal to step into a habitat with
# predator = 0 predators/100 km^2 than predator = 5 predators/100 km^2?
x1 <- data.frame(forage = 600, predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

x2 <- data.frame(forage = 600, predator = 5, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

# We can also ask for confidence intervals
logRSS <- log_rss(m1, x1 = x1, x2 = x2, ci = "se", ci_level = 0.95)

# RSS
exp(logRSS$df$log_rss)
# Confidence interval
exp(c(logRSS$df$lwr, logRSS$df$upr))

# Our animal is about 2.3x more likely to step into a habitat with
# predator = 0 than predator = 5. The 95% CI for that estimate is
# ~ 0.5 -- 11.0.

## Figure
# If we want to use this to make a figure, we can pass a sequence
# of values to x1. Remember, x2 must always be 1 row. Let's 
# visualize the RSS for forage vs. mean forage.
x1 <- data.frame(forage = seq(0, 800, length.out = 100), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

x2 <- data.frame(forage = mean(values(hab$forage)), 
                 predator = 0, 
                 sl_ = 50, log_sl_ = log(50),
                 cos_ta_ = 0)

logRSS <- log_rss(m1, x1, x2, ci = "se", ci_level = 0.95)

# We have a plot method for 'log_rss' objects to make a very basic figure.
plot(logRSS)

# But if we want more control, we can use ggplot with the 'df' data.frame.
ggplot(logRSS$df, aes(x = forage_x1, y = exp(log_rss), 
                      ymin = exp(lwr), ymax = exp(upr))) +
  geom_ribbon(color = "black", fill = "gray80", linetype = "dashed") +
  geom_line() +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  xlab("Forage at x1") +
  ylab("RSS vs. Mean Forage") +
  theme_bw()

# More complex iSSF ----
# (I.e., revisiting Step 3)
# Now let's try a much more complex model with interactions that add
# interesting biological realism.

# Suppose our animal is nocturnal. We know that it likes to sleep
# in forests during the day and forage in grasslands at night. So 
# we want selection for our cover types to vary with the time of day.

# Since our animal is nocturnal, we might also expect it's movement
# patterns to be different between day and night. So we want our
# movement parameters to interact with time of day.

# We can capture this complexity with interactions.

m2 <- covs2 %>% 
  # Make your own dummy variables -- R makes too many levels with
  # interactions between two categorical variables for clogit models.
  # Grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_end_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_end_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_end_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_end_ == "night")) %>%
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
summary(m2)

# ... visualize movement ----
# You cannot simply use 'update_sl_distr()' or 'update_ta_distr()' if
# you have interactions with your movement parameters. You need to pass 
# the fully predicted betas to the correct updating function, instead.

# Grab all the betas
b <- coef(m2)

# Update step-length distribution

# For tod = "day"

b_sl_day <- b[["tod_start_day:sl_"]]
b_log_sl_day <- b[["tod_start_day:log_sl_"]]

sl_distr_day <- update_gamma(sl_distr(m2),
                             beta_sl = b_sl_day,
                             beta_log_sl = b_log_sl_day)

# And what if tod = "night"?

b_sl_night <- b[["tod_start_night:sl_"]]
b_log_sl_night <- b[["tod_start_night:log_sl_"]]

sl_distr_night <- update_gamma(sl_distr(m2),
                               beta_sl = b_sl_night,
                               beta_log_sl = b_log_sl_night)

# Compile into a data.frame for easier plotting
sl_tod <- data.frame(tod = c("day", "night", "tent"),
                     shp = c(sl_distr_day$params$shape, 
                             sl_distr_night$params$shape,
                             sl_df$shp[1]),
                     scl = c(sl_distr_day$params$scale, 
                             sl_distr_night$params$scale,
                             sl_df$scl[1]))

# Plot
expand.grid(sl = seq(1, 1000, length.out = 100),
            tod = c("day", "night", "tent")) %>% 
  left_join(sl_tod) %>% 
  mutate(y = dgamma(sl, shape = shp, scale = scl)) %>% 
  ggplot(aes(x = sl, y = y, color = tod, linetype = tod)) +
  geom_line(size = 1) +
  scale_color_manual(name = "Time of Day", 
                     breaks = c("day", "night", "tent"),
                     labels = c("Day", "Night", "Tentative"),
                     values = c("goldenrod", "navy", "gray70")) +
  scale_linetype_manual(name = "Time of Day",
                        breaks = c("day", "night", "tent"),
                        labels = c("Day", "Night", "Tentative"),
                        values = c("solid", "solid", "dotdash")) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  theme_bw()

# As we hypothesized, our animal has longer steps during the night
# than during the day.

# Notice that both of these distributions are quite different from the 
# tentative distribution.

# Now how about the turn-angle distribution? 

cos_ta_day <- b[["cos_ta_"]]
ta_distr_day <- update_vonmises(ta_distr(m2), 
                                beta_cos_ta = cos_ta_day)

cos_ta_night <- b[["cos_ta_"]] + b[["cos_ta_:tod_start_night"]]
ta_distr_night <- update_vonmises(ta_distr(m2), 
                                  beta_cos_ta = cos_ta_night)

# Compile into a data.frame for easier plotting
ta_tod <- data.frame(tod = c("day", "night", "tent"),
                     k = c(ta_distr_day$params$kappa, 
                             ta_distr_night$params$kappa,
                             ta_df$k[1]))

# Plot
expand.grid(ta = seq(-pi, pi, length.out = 100),
            tod = c("day", "night", "tent")) %>% 
  left_join(ta_tod) %>% 
  rowwise() %>% 
  mutate(y = circular::dvonmises(ta, mu = 0, kappa = k)) %>% 
  ggplot(aes(x = ta, y = y, color = tod, linetype = tod)) +
  geom_line(size = 1) +
  scale_color_manual(name = "Time of Day", 
                     breaks = c("day", "night", "tent"),
                     labels = c("Day", "Night", "Tentative"),
                     values = c("goldenrod", "navy", "gray70")) +
  scale_linetype_manual(name = "Time of Day",
                        breaks = c("day", "night", "tent"),
                        labels = c("Day", "Night", "Tentative"),
                        values = c("solid", "solid", "dotdash")) +
  scale_x_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
  coord_cartesian(ylim = c(0, 0.4)) +
  xlab("Turn Angle (radians)") +
  ylab("Probability Density") +
  theme_bw()

# Notice, again, that both of these distributions are quite different from the 
# tentative distribution.

# ... visualize habitat selection ----
# This is a complex model with a lot of interesting relationships to
# visualize. Let's pick just one -- we'll visualize how selection for
# the different land cover types changes with time of day.

# Note that factors must have the same levels as in the fitted model.

# We want to compare each land cover category to "grassland". So the only
# thing that should change between x1 and x2 should be 'cover'. However,
# we want to do this for both day and night -- we need to make these
# comparisons separately, and then combine the results later for our
# figure.

# Day
x1_lc_day <- data.frame(forage = mean(values(hab$forage)),
                        temp = mean(values(hab$temp)),
                        predator = mean(values(hab$predator)),
                        cover = factor(c("grassland", "forest", "wetland"),
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("day", levels = c("day", "night")),
                        tod_end_ = factor("day", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))

x2_lc_day <- data.frame(forage = mean(values(hab$forage)),
                        temp = mean(values(hab$temp)),
                        predator = mean(values(hab$predator)),
                        cover = factor("grassland",
                                       levels = c("grassland", "forest", "wetland")),
                        tod_start_ = factor("day", levels = c("day", "night")),
                        tod_end_ = factor("day", levels = c("day", "night")),
                        sl_ = 100,
                        log_sl_ = log(100),
                        cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
# Calculate log-RSS
log_rss_lc_day <- log_rss(m2, x1 = x1_lc_day, x2 = x2_lc_day, ci = "se")

# Night
x1_lc_night <- data.frame(forage = mean(values(hab$forage)),
                          temp = mean(values(hab$temp)),
                          predator = mean(values(hab$predator)),
                          cover = factor(c("grassland", "forest", "wetland"),
                                         levels = c("grassland", "forest", "wetland")),
                          tod_start_ = factor("night", levels = c("day", "night")),
                          tod_end_ = factor("night", levels = c("day", "night")),
                          sl_ = 100,
                          log_sl_ = log(100),
                          cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
x2_lc_night <- data.frame(forage = mean(values(hab$forage)),
                          temp = mean(values(hab$temp)),
                          predator = mean(values(hab$predator)),
                          cover = factor("grassland",
                                         levels = c("grassland", "forest", "wetland")),
                          tod_start_ = factor("night", levels = c("day", "night")),
                          tod_end_ = factor("night", levels = c("day", "night")),
                          sl_ = 100,
                          log_sl_ = log(100),
                          cos_ta_ = 1) %>% 
  # Make your own dummy variables
  # Remember grassland is our reference level
  mutate(forest_day = as.numeric(cover == "forest" & tod_start_ == "day"),
         wetland_day = as.numeric(cover == "wetland" & tod_start_ == "day"),
         forest_night = as.numeric(cover == "forest" & tod_start_ == "night"),
         wetland_night = as.numeric(cover == "wetland" & tod_start_ == "night"))
# Calculate log-RSS
log_rss_lc_night <- log_rss(m2, x1 = x1_lc_night, x2 = x2_lc_night, ci = "se")

# Combine and plot
bind_rows(log_rss_lc_day$df, log_rss_lc_night$df) %>% 
  ggplot(aes(x = cover_x1, y = log_rss, color = tod_end_x1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(0.25)) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), 
                position = position_dodge(0.25), width = 0.1) +
  scale_color_manual(name = "Time of Day",
                     breaks = c("day", "night"),
                     values = c("goldenrod", "navy")) +
  xlab("Land Cover") +
  ylab("log-RSS vs Grassland") +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw()

# We can see that, relative to grassland, forest is selected during the day
# and avoided at night. Wetland is avoided during both day and night.
