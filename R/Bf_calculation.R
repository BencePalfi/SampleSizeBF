#'BF function based on Dienes & McLatchie, 2018
#'
#'Modified so that H1 is represented by normal distribution (rather than t), hence there is no 'dftheory' argument
#'
#'@param sd SE of data
#'
#'@section notes: this function is based on the Dienes & Mclatchie (2018) calculator
#'
#'@export
Bf_normal <- function(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail = 2) { 
  area <- 0
  normarea <- 0
  theta <- mean_of_theory - 5 * sd_of_theory
  incr <- sd_of_theory / 200
  for (A in -1000:1000) {
    theta <- theta + incr
    dist_theta <- dnorm((theta - mean_of_theory) / sd_of_theory)
    if (identical(tail, 1)) {
      if (theta <= mean_of_theory) {
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained - theta) / sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta * incr
  }
  LikelihoodTheory <- area / normarea
  Likelihoodnull <- dt(obtained / sd, df = dfdata)
  BayesFactor <- LikelihoodTheory / Likelihoodnull
  BayesFactor
}

#'@export
Bf_uniform <- function(sd, obtained, dfdata, lower, upper) {
  area <- 0
  normarea <- 0
  theta <- lower
  range <- upper - lower
  incr <- range / 2000
  for (A in -1000:1000) {
    theta <- theta + incr
    dist_theta <- 1 / range
    height <- dist_theta * dt((obtained - theta) / sd, df = dfdata)
    area <- area + height * incr
  }
  LikelihoodTheory <- area
  Likelihoodnull <- dt(obtained / sd, df = dfdata)
  LikelihoodTheory / Likelihoodnull
}

#'@export
Bf_cauchy <- function(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail = 2) {
  area <- 0
  normarea <- 0
  theta <- mean_of_theory - 5 * sd_of_theory
  incr <- sd_of_theory / 200
  for (A in -1000:1000) {
    theta <- theta + incr
    dist_theta <- dt((theta - mean_of_theory) / sd_of_theory, df = 1)
    if(identical(tail, 1)) {
      if (theta <= mean_of_theory) {
        dist_theta <- 0
      } else {
        dist_theta <- dist_theta * 2
      }
    }
    height <- dist_theta * dt((obtained-theta) / sd, df = dfdata)
    area <- area + height * incr
    normarea <- normarea + dist_theta * incr
  }
  LikelihoodTheory <- area / normarea
  Likelihoodnull <- dt(obtained / sd, df = dfdata)
  LikelihoodTheory/Likelihoodnull
}