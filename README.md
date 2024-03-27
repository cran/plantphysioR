# plantphysioR <img src="man/figures/logo.png" align="right" height="120"/>

The goal of plantphysioR is to provide basic functions related to plant physiology, it contains functions related to various yield related indices, pigment content calculations and PEG6000 calculations

## Installation

You can install the development version of plantphysioR from [GitHub](https://github.com/rameshram96/plantphysioR) with:

``` r
# install.packages("devtools")
devtools::install_github("rameshram96/plantphysioR")
```

## Example

``` r
library(plantphysioR)
#all_indices 
Ms<-mean(yield_data$Ys)
Ms<-mean(yield_data$Ys)
Yp<-yield_data$Yp
Ys<-yield_data$Ys
all_indices(Yp,Ys,Mp,Ms)
#This function calculates all the yield related indices inclued in the package
# calculation of Stress susceptibility index
ss_index(500, 350, 450, 370)
```

<!-- badges: start -->

[![R-CMD-check](https://github.com/rameshram96/plantphysioR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rameshram96/plantphysioR/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->
