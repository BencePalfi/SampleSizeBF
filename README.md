
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SampleSizeBF

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The goal of SampleSizeBF is TODO: BENCE

## Features

-   Calculating Bayes factor with normal, Cauchy, and uniform
    distributions as a prior
-   Samplesize planning for Bayes factor using the same priors
-   Simulate samples to verify the results of the methods presented in
    the package

## Usage

SampleSizeBF can be used either via the Shiny app or via R.

## Using the Shiny app

You can use the app at [LINK](http://www.staggeringbeauty.com/)

You can alternatively run the app locally on your own computer by
following these instructions:

Install the development version (SampleSizeBF is not available from
CRAN) from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("BencePalfi/SampleSizeBF")
```

Running the app.

``` r
SampleSizeBF::run_app()
```

You can read more on how to use the SampleSizeBF app in
vignette(“app\_use”).

## Installation

You can install the released version of SampleSizeBF from with:

``` r
install.packages("SampleSizeBF")
```

## Using the package

You can install the development version of the package from GitHub by
running the code presented in the previous section.

You can read more on how to use the SampleSizeBF package to create
reports from R in vignette(“local\_use”).
