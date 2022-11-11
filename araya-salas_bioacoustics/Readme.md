# Bioacoustic data analysis

### November 7, 2022

## [Marcelo Araya-Salas](marceloarayasalas.weebly.com/)
#### Neuroscience Research Center/Biology Department, University of Costa Rica

---

The files shared here aim to provide an overview of the tools available in R for streamlining bioacoustic data analysis. 

---

## R Packages

The webinar focuses on the R packages [warbleR](https://marce10.github.io/warbleR/), [ohun](https://marce10.github.io/ohun/), [Rraven](https://marce10.github.io/Rraven/) and [PhenotypeSpace](https://marce10.github.io/PhenotypeSpace/) for preparing our acoustic data and quantifying the structure of sounds. 

Note that the packages [ohun](https://marce10.github.io/ohun/) and [PhenotypeSpace](https://marce10.github.io/PhenotypeSpace/) haven't been submitted to CRAN yet (nov 7 2022) and must be installed from github.

The latest version of the packages needed to run the code can be installed as follows:

``` r
remotes::install_github("maRce10/Rraven")
remotes::install_github("maRce10/warbleR")
remotes::install_github("maRce10/ohun")
remotes::install_github("maRce10/PhenotypeSpace")

```

---

### Scripts
  
The R script can be downloaded [from here](https://raw.githubusercontent.com/eco4cast/Statistical-Methods-Seminar-Series/main/araya-salas_bioacoustics/example_code_bioacoustics_Araya-Salas_2022.Rmd) and its also found within the subfolder "araya-salas_bioacoustics" 
  
### Data

All data used in the example R code is found in the "data" subfolder. Download that folder into the same directory in which the Rmd script is found to run the code. 

### Webinar

The webinar was recorded and is available on youtube [in this link](https://youtu.be/dYlkDCUTbAs).
