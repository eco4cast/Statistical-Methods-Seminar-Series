# Spatial phylogenetics with {phylospatial} <a href="https://matthewkling.github.io/phylospatial/"><img src="images/phylospatial-logo.png" align="right" height="139"/></a>

Matthew Kling · ESA Statistical Ecology / EFI Statistical Methods Seminar Series · November 2, 2026

An introduction to spatial phylogenetics (mapping the distribution of evolutionary history) and to the {phylospatial} R package, using a dataset of 443 California moss species. The seminar covers phylogenetic diversity and endemism, null-model significance tests and CANAPE, phylogenetic turnover and regionalization, and conservation prioritization.

## Materials

- `phylospatial-webinar.qmd`: the slides as a Quarto file. Every slide's code runs, so you can follow along chunk by chunk.
- `phylospatial-webinar.html`: the rendered slides; download it and open in a browser.
- `data/moss_points.csv`: moss occurrence records used in the gridding example.
- `images/`: figures used in the slides.

## Following along

Install the packages used in the workshop:

``` r
install.packages(c("phylospatial", "nullcat", "ggplot2", "tidyterra", "patchwork",
                   "prioritizr", "highs"))
```

One section uses a helper function that is currently only in the development version of phylospatial:

``` r
install.packages("pak")
pak::pak("matthewkling/phylospatial")
```

Then open `phylospatial-webinar.qmd` and run chunks interactively, with this folder as your working directory (it reads `data/moss_points.csv` by relative path). 

## Links

- Package documentation and vignettes: <https://matthewkling.github.io/phylospatial/>
- Source code and issues: <https://github.com/matthewkling/phylospatial>
- Kling (2025) *Methods in Ecology and Evolution*: <https://doi.org/10.1111/2041-210X.70056>
- Contact: <matthew.kling@umontana.edu>
