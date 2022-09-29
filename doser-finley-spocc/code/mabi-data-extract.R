# Fit a trend model with EAWP for this data set, then use BBS for the multi-species stuff. 
rm(list = ls())
library(tidyverse)
library(sf)

load("data/mr-ms-rs-data.R")
# Eastern Wood Pewee
sp.indx <- 10
y.small <- y[2, sp.indx, 1:J[2], , ]
regen.small <- regen[2, 1:J[2]]
basalArea.small <- basalArea[2, 1:J[2]]
percentForest.small <- percentForest[2, 1:J[2]]

# Reformat the y data
for (i in 1:nrow(y.small)) {
  print(i)
  for (j in 1:ncol(y.small)) {
    for (k in 1:dim(y.small)[3]) {
      if (!is.na(y.small[i, j, k])) {
        if (y.small[i, j, k] > 0 & k < dim(y.small)[3]) {
          y.small[i, j, (k + 1):dim(y.small)[3]] <- NA
        }
      }
      y.small[i, j, k] <- ifelse(y.small[i, j, k] > 0, 1, 0) 
    } # k
  } # j
} # i
# Read in coordinates from DFWZ20
coord.df <- read.csv("~/Dropbox/DFWZ20/data/mabi-map/Points.csv")
coords <- coord.df %>%
  filter(Admin_Unit_Code == 'MABI')  %>%
  select(Longitude, Latitude)

coords.sf <- st_as_sf(coords, 
		      coords = c('Longitude', 'Latitude'), 
		      crs = '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
coords.sf <- coords.sf %>%
  st_transform(crs = "+proj=utm +zone=18 +datum=NAD83 +units=m +no_defs")
coords.proj <- st_coordinates(coords.sf)

year.cov <- matrix(year.regions[2, ], nrow(coords.sf), ncol(year.regions), byrow = TRUE)

occ.covs <- list(regen = regen.small, 
		 basalArea = basalArea.small, 
		 percentForest = percentForest.small, 
		 site.index = 1:nrow(coords.sf), 
                 year = year.cov)
det.covs <- list(basalArea = basalArea.small, year = year.cov)
data.list <- list(y = y.small, occ.covs = occ.covs, det.covs = det.covs,
		  coords = coords.proj)
attr(data.list$occ.covs$regen, 'names') <- NULL
attr(data.list$occ.covs$basalArea, 'names') <- NULL
attr(data.list$occ.covs$percentForest, 'names') <- NULL
attr(data.list$det.covs$basalArea, 'names') <- NULL

save(data.list, file = 'data/doser2021EcoApps.rda')
