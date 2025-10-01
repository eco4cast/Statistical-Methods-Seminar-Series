# The `unmarked` R package for fitting hierarchical models of animal abundance and occurrence

### [Ken Kellner](https://kenkellner.com)

#### Michigan State University

## About

Imperfect detection is a common problem in ecological sampling which can lead to biased estimates of site occupancy and abundance.
`unmarked` is an R package for fitting occupancy and abundance models that account for imperfect detection.
The package provides a workflow for fitting models to a wide variety of data types and experimental designs. 

## Talk outline

1. Overview of hierarchical models
2. The `unmarked` package
    - Organization
    - Supported models
3. Package demonstration
4. Recent features
    - Random effects
    - Simulation
    - Power analysis
5. Future directions
6. Questions

## Installing `unmarked`

Install from CRAN with:

```r
install.packages("unmarked")
```

The example code also uses the `terra` library:

```r
install.packages("terra")
```

## Repository structure

`unmarked_presentation.ipynb`: The presentation as a Jupyter notebook

`unmarked_presentation.R`: All R code used in the presentation

`unmarked_presentation.Rmd`: Jupyter notebook converted to an Rmarkdown file

`unmarked_presentation.html`: Presentation as html

`Makefile`: Describes how the files above were derived from the Jupyter notebook

### data

`y.csv`: Detection-nondetection data for crossbill

`site_covs.csv`: Site-level covariates (elevation, forest)

`dates.csv`: Sampling (ordinal) dates corresponding to each element in `y`

`example_data.R`: Script to create the CSVs from the original `crossbill` dataset provided in `unmarked`

### images

Image files used in the presentation

## More resources

https://cran.r-project.org/package=unmarked (many vignettes)

https://github.com/ecoverseR/unmarked

Previous ESA/EFI seminar on [multispecies occupancy models in unmarked](https://htmlpreview.github.io/?https://github.com/eco4cast/Statistical-Methods-Seminar-Series/blob/main/Rota-MSOM/MSOM-Presentation.html)
