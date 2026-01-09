# Building Spatial Statistical Models in **R** Using `spmodel` <a href="https://usepa.github.io/spmodel/"><img src="https://github.com/USEPA/spmodel/blob/main/man/figures/logo.png" align="right" height="240" alt="spmodel website" /></a> 

### Michael Dumelle

#### United States Environmental Protection Agency

## Repository Contents

* dumelle-spmodel.html: The presentations slides (reveal.js format).
* dumelle-spmodel.qmd: The Quarto document generating dumelle-spmodel.html.
* references.bib: The BibTeX references in dumelle-spmodel.html.
* figures: The figures in dumelle-spmodel.html.
* _extensions: The Quarto extension for the slide pointer in dumelle-spmodel.html (by pressing the "q" key).
* dumelle-spmodel.R: The R code in dumelle-spmodel.qmd (the output of `knitr::purl("dumelle-spmodel.qmd")`).
* README.md: This README.

## **R** Packages Used

The following **R** packages are required to compile `dumelle-spmodel.qmd` and can be installed by running:

```{r}
install.packages("tidyverse") # for general data wrangling
install.packages("broom") # for tidy() of lm() objects
install.packages("kableExtra") # for tables
install.packages("sf") # for spatial data operations
install.packages("spmodel") # for spatial modeling
install.packages("tigris") # for southwestern state boundaries
```

## Abstract

Statistical and machine learning models often assume observations are independent of one another. This assumption is typically impractical for spatial data, as nearby observations tend to be more similar than distant observations. Models that ignore spatial dependence can suggest misleading relationships between response and explanatory variables, which obfuscates ecological understanding and leads to ineffective management strategies. The `spmodel` (https://usepa.github.io/spmodel/) R package is new, freely available software for building spatial models. `spmodel`’s `splm()` and `spglm()` functions are used to fit linear and generalized linear models which formally account for spatial dependence among nearby observations. The `splm()` and `spglm()` functions build upon R’s familiar `lm()` and `glm()` functions, making the transition from nonspatial models to spatial models relatively seamless. In this talk, Dr. Dumelle will provide an overview of spatial statistics and show how to use `spmodel` to fit, interpret, and make predictions using spatial models.

Dr. Michael Dumelle is the lead statistician for the United States (US) Environmental Protection Agency’s National Aquatic Resource Surveys (NARS), a national monitoring program designed to characterize the physical, chemical, and biological condition of aquatic ecosystems in the US. Dr. Dumelle’s research interests are in spatial statistics, missing data, spatial survey design, and software development. He is the maintainer of the R packages `spmodel` for building spatial statistical models, `SSN2` for building spatial statistical models specifically on stream networks, and `spsurvey` for spatially balanced random sampling.

