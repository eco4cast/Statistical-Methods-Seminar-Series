################################################################X
#-------------ESA Statistical Methods Seminar Series------------X
#------Avgar & Smith -- Integrated Step Selection Analysis------X
#-------------------------07 March 2022-------------------------X
################################################################X

# Extra script: simulating habitat rasters
# Brian J. Smith <brian.smith@usu.edu>

# Load packages----
library(tidyverse)
library(raster)
library(NLMR)

# Generate habitats ----

# ... mean annual temp ----
# units are degrees C
# reasonable values are 0 - 20 °C

h1 <- nlm_fbm(ncol = 300, nrow = 300, 
              resolution = 50, fract_dim = 1.2,
              user_seed = 20220203 + 1)

# Stretch to range
temp <- (values(h1)) * 20

# Check
hist(temp, breaks = 30)

# Place in raster
temp_rast <- setValues(h1, temp)

# ... forage ----
# units are g/m2
# reasonable values are 0 - 1000
# let's make it 0 inflated (rocky outcrops?)

h2 <- nlm_fbm(ncol = 300, nrow = 300, 
              resolution = 50, fract_dim = 0.65,
              user_seed = 20220203 + 2)

# Stretch to range (create some negatives)
forage <- (values(h2) - 0.2)/0.8  * 1000

# Cut any negative values to 0
forage[forage < 0] <- 0

# Check
hist(forage, breaks = 30)

# Place in raster
forage_rast <- setValues(h1, forage)

# ... predation risk ----
# units are predators/100 km2
# reasonable values are 0 - 15 (think large predator like wolf)
# want mean to vary through space

# Centroids
cents <- lapply(1:8, function(x) {
  set.seed(x * 3)
  x <- runif(1, 0, 15000)
  y <- runif(1, 0, 15000)
  return(data.frame(x = x, y = y))
}) %>% 
  bind_rows() %>% 
  as.matrix()

h3cent <- distanceFromPoints(h1, cents)
h3cent <- h3cent/max(values(h3cent)) * 2

h3grad <- nlm_distancegradient(ncol = 300, nrow = 300, resolution = 50,
                            origin = c(120, 170, 10, 10))

# Combine
h3 <- sum(h3cent, h3grad)

# Scale
h3 <- max(values(h3)) - h3

# Place in raster
pred_rast <- h3/max(values(h3)) * 15

# ... landcover ----
# categories are:
#   - grassland (50%)
#   - forest (30%)
#   - wetland (20%)

set.seed(20220105 + 4)
h4 <- nlm_mosaictess(ncol = 300, nrow = 300, 
                     resolution = 1, germs = 400)

germ2lc <- data.frame(value = sort(unique(values(h4)))) %>% 
  mutate(landcover = case_when(
    value < 0.5 ~ 1,
    value < 0.8 ~ 2,
    TRUE ~ 3
  ))

lc_vals <- data.frame(value = values(h4)) %>% 
  left_join(germ2lc)

lc <- setValues(h1, lc_vals$landcover)

# Check
table(values(lc))/(300*300)

# Stack and convert to UTM ----
rast <- stack(forage_rast, temp_rast, pred_rast, lc)
names(rast) <- c("forage", "temp", "pred", "lc")
# Check
plot(rast)

# Get data and coordinates
hab_dat <- as.data.frame(rast, xy = TRUE)

# Shift coordinates
hab_dat$x <- hab_dat$x + 447000
hab_dat$y <- hab_dat$y + 4626000

# Back to raster
hab <- rasterFromXYZ(hab_dat, res = 50, crs = 32612)

# Save ----
dir.create("data", showWarnings = FALSE)
writeRaster(hab, "data/habitat.tif", overwrite = TRUE)
