rm(list = ls())
library(spOccupancy)
library(tidyverse)
library(sf)

load("~/Dropbox/DKFSWZ22/data/spOcc-bbs-data.rda")


# Species info
sp.names <- dimnames(data.list$y)[[1]]

y.sums <- sort(apply(data.list$y[, , 20, ], 1, sum, na.rm = TRUE))

my.sp <- c('OVEN', 'SCTA', 'AMRE', 'WBNU', 'VEER', 'BAWW', 'PIWO', 
	   'BHVI', 'PIWA', 'BTNW', 'HAWA', 'CORA')
sp.indx <- which(sp.names %in% my.sp)

y <- data.list$y[sp.indx, , 20, ]

site.indx <- which(apply(y, 2, sum, na.rm = TRUE) > 0)
y <- y[, site.indx, ]
occ.covs <- data.frame(elev = data.list$occ.covs$elev[site.indx], 
		       bcr = data.list$occ.covs$BCR[site.indx])
det.covs <- list(day = data.list$det.covs$day[site.indx, 20], 
		 obs = data.list$det.covs$obs[site.indx, 20])
coords <- data.list$coords[site.indx, ]

data.bbs <- list(y = y, 
		 occ.covs = occ.covs, 
		 det.covs = det.covs,
		 coords = coords)

save(data.bbs, file = "data/bbs2019.rda")

# Single species data -----------------
# Veery
y.one <- y[which(my.sp == 'VEER'), , ]
data.veery <- list(y = y.one, 
		   occ.covs = occ.covs, 
		   det.covs = det.covs, 
		   coords = coords)
save(data.veery, file = 'data/bbs2019Veery.rda')

# Get prediction data -----------------------------------------------------
load("~/Dropbox/ZL22/ICOM/data/ne-data-bundle.rda")
elev.pred <- occ.covs$elev
elev.pred <- ifelse(is.na(elev.pred), 0, elev.pred)
centers <- st_centroid(grid.ne)
coords.0 <- st_coordinates(centers) / 1000
save(elev.pred, coords.0, file = 'data/bbs-2019-pred-data.rda')

