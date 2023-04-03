MixSIAR tutorial
=============

[MixSIAR](http://brianstock.github.io/MixSIAR/index.html) is an R package that helps you create and run Bayesian mixing models to analyze biotracer data (i.e. stable isotopes, fatty acids), following the [MixSIAR model framework](https://peerj.com/articles/5096/).

## Install

1. Download and install/update [R](https://cran.r-project.org/).

2. Download and install [JAGS](http://mcmc-jags.sourceforge.net/).

3. Install MixSIAR. Open R and run:
```
install.packages("MixSIAR", dependencies=TRUE)
library(MixSIAR)
```

I recommend installing the GitHub version because it has the latest changes and bug fixes not yet on CRAN. This is the version I will use today.
```
remotes::install_github("brianstock/MixSIAR", dependencies=T)
```

## Ex 1: Wolves

http://brianstock.github.io/MixSIAR/articles/wolves_ex.html

## Ex 2: Modifying MixSIAR output

http://brianstock.github.io/MixSIAR/articles/modify_output.html

## Full tutorial

Several more [vignettes](http://brianstock.github.io/MixSIAR/articles/index.html) explore additional features of MixSIAR.

There is also an extensive user manual included in the package install. To find the directory location on your computer:
```
find.package("MixSIAR")
```

Alternatively, you can download the manual from the GitHub site [here](https://github.com/brianstock/MixSIAR/blob/master/inst/mixsiar_manual_small.pdf).

Clean, runnable `.R` scripts for each vignette are also available in the `example_scripts` folder of the `MixSIAR` package install:
```
library(MixSIAR)
mixsiar.dir <- find.package("MixSIAR")
file.path(mixsiar.dir, "example_scripts")
```

You can then run the Wolves example script with:
```
setwd("choose/where/to/save/output")
source(file.path(mixsiar.dir, "example_scripts", "mixsiar_script_wolves.R"))
```
