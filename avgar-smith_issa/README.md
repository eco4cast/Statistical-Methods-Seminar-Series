# Integrated Step Selection Analysis (iSSA)

### March 7, 2022

## Tal Avgar and Brian J. Smith

### Department of Wildland Resources, Utah State University

---

## About iSSA

A habitat selection function is a model of the relative probability that an available spatial unit will be used by an animal given its habitat value, but how do we appropriately define availability? In an integrated Step-Selection Analysis (iSSA), availability is defined by the animal’s ‘selection-free movement kernel’, which is fitted in conjunction with a conditional habitat-selection function. Parameter estimates are obtained using a conditional-logistic regression by contrasting each ‘used step’ (a straight line connecting two consecutive observed positions of the animal) against a set of ‘available steps’ (randomly sampled from one of several possible theoretical distributions). iSSA thus relaxes the implicit assumption that movement is independent of habitat selection and instead allows simultaneous inference on both processes, resulting in an empirically parameterized mechanistic space-use model.

---

## R Packages

We will be focusing on the R package `amt` for preparing our data, running an iSSA, and interpreting the results.

The latest version of `amt` is available on CRAN. [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/amt)](https://cran.r-project.org/package=amt) 

You can find the development version on [GitHub](https://github.com/jmsigner/amt). You will need additional R build tools if you want to install the development version (i.e., Rtools on Windows; Xcode on Mac).

You might run into some warnings (or worse) if you have older versions of other packages. In particular, it would help to update `sf` and `raster` before working through this code.

---

## Repository Structure

### Presentation

A PDF of the slides from the main presentation is saved as `iSSA_ESAwebinar.pdf`. The FAQ slides are at the end of the main slides.

### Q & A

A follow-up document with more questions and answers is available as `Q_and_A.md`. It should be rendered properly here on GitHub.

### Scripts

This repository contains 4 R scripts. The main R script, `amt_demo.R`, is what we work through during the webinar. In addition, there are supplementary scripts contained in the directory `extra/`:

  1. `sim_habitat.R` -- simulates the habitat rasters used for the analysis
  2. `sim_tracks.R` -- simulates the animal location data from an iSSA
  3. `quick_issa.R` -- fits the full iSSF from main demo in piped workflow
  
### Data

The data simulated in extra scripts 1 & 2 are saved in the directory `data/`:

  1. `data/habitat.tif` -- 4-band GeoTIFF containing habitat variables: 
  
    a. Forage (g/m^2^) -- treated as a *resource* (more is better)
    b. Mean temperature (°C) -- treated as a *condition* (intermediate value is optimal)
    c. Predator density (predators/100 km^2^) -- treated as a *risk* (less is better)
    d. Landcover (categories) -- grassland, forest, or wetland; no *a priori* designation as resource, condition, or risk. 
    
  2. `data/gps.csv` -- comma delimited spreadsheet with 4 columns:
  
    a. ID -- Unique identifier for the individual animal (only 1 individual)
    b. utm_e -- Easting in UTM Zone 12 (EPSG: 32612)
    c. utm_e -- Northing in UTM Zone 12 (EPSG: 32612)
    d. timestamp -- Date and time in timezone "US/Mountain"
    
---

## Updates to `amt` Coming Soon
We didn't cover these topics in the webinar, partly for time, and partly because we don't have polished functions in `amt` for them just yet. They are currently in the works, and you can expect to see them in the coming weeks and months.

  * Used-habitat calibration plots (UHC plots; [Fieberg et al. 2018](https://doi.org/10.1111/ecog.03123))
  * Simulators for fitted models (e.g., [Signer et al. 2017](https://doi.org/10.1002/ecs2.1771))
  
---
