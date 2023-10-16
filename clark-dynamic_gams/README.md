# Dynamic Generalized Additive Models with `mvgam`

### Nicholas J Clark

## About
Time series analysis and forecasting are standard goals in applied ecology. But ecological forecasting is difficult because ecology is complex. The abundances of species, for example, fluctuate for many reasons. Food and shelter availability limit survival. Biotic interactions affect colonization and vital rates. Severe weather events and climate variation alter habitat suitability. These sources of variation make it difficult to understand, let alone predict, ecosystem change. Moreover, most available time series software cannot handle features that dominate ecological data, including overdispersion, clustering, missingness, discreteness and nonlinear effects.

In this talk, I will introduce Dynamic Generalized Additive Models (DGAMs) as one solution to meet this complexity. I will illustrate a number of models that can be tackled with the [{mvgam} R package](https://nicholasjclark.github.io/mvgam/), which builds Stan code to specify probabilistic Bayesian models that include nonlinear smooth functions, random effects and dynamic processes, all with a simple interface that is familiar to most R users.

## Slides
The slidedeck can be [accessed as an html version](https://nicholasjclark.github.io/EFI_seminar/EFI_talk_slidedeck#1) or [downloaded as a PDF](https://github.com/nicholasjclark/EFI_seminar/raw/main/EFI_talk_slidedeck.pdf)

## System requirements

Please be sure to have at least version 4.1 &mdash; *and preferably version 4.2* &mdash; of `R` installed. We will make use of several `R` packages that you'll need to have installed. Prior to the start of the course, please run the following code to update your installed packages and then install the required packages:

```r
# update any installed R packages
update.packages(ask = FALSE, checkBuilt = TRUE)

# packages to install
pkgs <- c("dplyr", "gratia", "ggplot2",
          "marginaleffects", "tidybayes", "zoo",
          "viridis", "remotes", "brms")

# install those packages by setting Ncpus to number of CPU cores you have available
install.packages(pkgs, Ncpus = 2)
```

For fitting models in `mvgam`, *it is highly recommended that you use the `Cmdstan` backend*, with the `cmdstanr` interface, rather than using `rstan`. To install `Cmdstan` and the relevant interface, first install `cmdstanr` using the following:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
```

And then [follow instructions provided by the `Stan` development team here](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) to ensure the backend is installed and the toolchain is setup properly


Finally, we will make use of the development version of the `mvgam` package as it is not quite ready for CRAN. You can install this package using the following:

```r
# Download and install mvgam
remotes::install_github('nicholasjclark/mvgam', force = TRUE)
```
