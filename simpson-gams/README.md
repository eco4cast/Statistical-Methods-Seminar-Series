# Generalized Additive Models (GAMs)

## Gavin Simpson

### Department of Animal Science, Aarhus University, Denmark

## About

Generalized Additive Models were introduced as an extension to linear and generalized linear models, where the relationships between the response and covariates are not specified up-front by the analyst but are learned from the data themselves. This learning is achieved by viewing the effect of a covariate on the response as a smooth function, rather than following a fixed form (linear, quadratic, etc). The smooth functions are represented in the GAM using penalized splines, in which a penalty against fitting overly-complex functions is employed. GAMs are most useful when the relationships between covariates and response are non linear, and GAMs have found particular use for modelling *inter alia* spatiotemporal data.

The presentation will briefly explain what a GAM is and how penalized splines work before focusing on the practical aspects of fitting GAMs to data using the [mgcv R package](https://cran.r-project.org/package=mgcv), and will be most useful to ecologists who already have some familiarity with linear and generalized linear models.

## R Packages

We need to use the development version of {gratia} for some aspects of today's seminar. There has been an annoying bug in `gratia::rootgram()` in the development version of {gratia} on GitHub for a while, however. Even if you used {gratia} previously, you'll need to reinstall from GitHUb

```r
# install.packages("remotes", Ncpus = 2)
remotes::install_github("gavinsimpson/gratia")
```

The following packages are also needed for the webinar; many will be installed already if you have installed {gratia}. Run the following code to check what is installed. Update packages if you have difficulties with the tidyverse:

```r
# change Ncpus to match your available CPU cores
update.packages(checkBuilt = TRUE, ask = FALSE, Ncpus = 2)
```
```r
pkgs <- c("here", "readr", "janitor", "mgcv", "dplyr",
          "ggplot2", "ggrepel", "mvnfast", "patchwork", "tidyr",
          "ppgam")
find.package(pkgs)
```

Check that there are no warnings. If there are warnings, install any packages listed in the output from `find.package()`. **Note**, don't reinstall {gratia} at this when you do this as we need the development version from GitHub.