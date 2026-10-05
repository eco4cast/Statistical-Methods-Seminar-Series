# Spatial Capture-Recapture Intro



This webinar will cover a basic introduction to spatial capture-recapture methods.  We will implement the models in R using the NIMBLE package, so a working knowledge of NIMBLE (or jags/bugs) is assumed but we will cover the basics. 



To follow along with the coding demonstration, you will need the R package NIMBLE. Please review the instructions here, https://r-nimble.org/download.html. You will also need a working C++ complier and the igraph, coda, R6, pracma, and numDeriv packages, see section 4 of the NIMBLE user manual for details: https://r-nimble.org/manual/cha-installing-nimble.html  



We will also use the coda and MCMCvis packages in this demonstration, and for the habitat script, you will need the terra, raster, and reshape2 packages.
 - install.packages(c("coda", "MCMCvis", "terra", "raster", "reshape2"))

**Still under construction**

**Data** folder contains the data needed for all of the exercises: 

- beardata.rda which has the following:

  - trapmat - trap id, x location, y location (UTM coordinates divided by 1000, units are kms)

  - bearArray - 3D array of captures by individual, trap, occasion

  - flat - encounter data formatted for the R package secr

  - sex – the biological sex assigned to individuals 

- forest_cover.tif
  - Percent forest cover from NLCD dataset for 2006. Data originally at 30m resolution, resampled here to 500m resolution. UTM Zone 18N. 

**Code** folder contains all the R scripts to run the exercises presented:

- ModelSCR\_0.R - fits an SCR model with no covariates to the bear data

- ModelSCR\_sex.R - fits an SCR model with sex as covariates on baseline detection and sigma to the bear data

- ModelSCR\_behavior.R - fits an SCR model with a one time behavior response to capture to the bear data
  
- ModelSCR\_habitat.R - fits an SCR model with an inhomogeneous point process where density varies as a function of forest cover



**Secr** folder (OPTIONAL): The secr folder contains script files to run all the models discussed in the webinar in the R package secr. I'm not going to cover the R package secr in this webinar; however, the secr package created and maintained by Murry Efford is a powerful and flexible tool for analyzing spatial capture recapture models with maximum likelihood.

