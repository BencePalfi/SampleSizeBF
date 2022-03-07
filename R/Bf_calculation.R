#' Calculate Bayes factor with normal prior
#' 
#' This function calculates the Bayes factor, i.e. the evidence for H1 over H0.
#' The function models the predictions of H1 with a normal distribution, and
#' assumes that the data generation process (i.e. the likelihood function) follows
#' a t-distribution. The function is based on the Dienes & Mclatchie (2018)
#' calculator. 
#' 
#' @family Bayes factor calculation functions
#'
#' @param sd numeric. Calculated SE of the observed effect. See \code{\link{calculate_se}}.
#' @param obtained numeric. The observed effect in raw units.
#' @param dfdata integer. The degrees of freedom of the sample. See \code{\link{calculate_df}}.
#' @param mean_of_theory numeric. The mode of the model of H1 (i.e. prior distribution).
#' @param sd_of_theory numeric. The standard deviation of the model of H1 (i.e. prior distribution).
#' @param tail integer. Either `1` for a half-normal prior distribution, or `2` for a two-tailed prior distribution.
#' 
#' @section Note: Please enter a positive obtained value if the observed effect
#'   is in the same directions as the effect predicted by H1, and a negative obtained value
#'   if not. This is critical to get a correct Bayes factor when you have a directional
#'   hypothesis (`tail = 1`). Note that the `sd_of_theory` must be a positive value.
#' 
#' @return The function returns the Bayes factor (evidence for H1 over H0) as a single
#' numeric value.
#' 
#' @examples 
#' # One-tailed test with the observed effect in the same direction as H1 predicted
#' SampleSizeBf::Bf_normal(sd = 0.4, obtained = 1, dfdata = 98, mean_of_theory = 0,
#'                         sd_of_theory = 1.5, tail = 1)
#'                         
#' # One-tailed test with the observed effect in the opposite direction as H1 predicted
#' SampleSizeBf::Bf_normal(sd = 0.4, obtained = -0.3, dfdata = 98, mean_of_theory = 0,
#'                         sd_of_theory = 1.5, tail = 1)
#' @export
Bf_normal <- function(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail = c(1, 2)) { 
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

#' Calculate Bayes factor with uniform prior
#'
#' This function calculates the Bayes factor, i.e. the evidence for H1 over H0.
#' The function models the predictions of H1 with a uniform distribution, and
#' assumes that the data generation process (i.e. the likelihood function) follows
#' a t-distribution. The function is based on the Dienes (2008) calculator. 
#'
#' @param sd numeric. Calculated SE of the observed effect. See \code{\link{calculate_se}}.
#' @param obtained numeric. The observed effect in raw units.
#' @param dfdata integer. The degrees of freedom of the sample. See \code{\link{calculate_df}}.
#' @param lower numeric. The smallest predicted effect size by H1.
#' @param upper numeric. The largest predicted effect size by H1.
#'
#' @section Note: Please set either `lower` or `upper` equal to 0 in order to
#'   model a directional H1.
#' 
#' @return The function returns the Bayes factor (evidence for H1 over H0) as a single
#' numeric value.
#' 
#' @examples 
#' # One-tailed test with the observed effect in the same direction as H1 predicted
#' SampleSizeBf::Bf_uniform(sd = 0.4, obtained = 1, dfdata = 98, lower = 0,
#'                         upper = 1.5, tail = 1)
#'                         
#' # One-tailed test with the observed effect in the opposite direction as H1 predicted
#' SampleSizeBf::Bf_uniform(sd = 0.4, obtained = -0.3, dfdata = 98, lower = 0,
#'                         upper = 1.5, tail = 1)
#'
#' @export
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

#' Calculate Bayes factor with Cauchy prior
#'
#' This function calculates the Bayes factor, i.e. the evidence for H1 over H0.
#' The function models the predictions of H1 with a Cauchy distribution, and
#' assumes that the data generation process (i.e. the likelihood function) follows
#' a t-distribution. The function is based on the Dienes & Mclatchie (2018)
#' calculator. 
#'
#' @param sd numeric. Calculated SE of the observed effect.
#' @param obtained numeric. The observed effect in raw units.
#' @param dfdata integer. The degrees of freedom of the sample.
#' @param mean_of_theory numeric. The mode of the model of H1 (i.e. prior distribution).
#' @param sd_of_theory numeric. The scale parameter of the model of H1 (i.e. prior distribution).
#' @param tail integer. Either `1` for a half-normal prior distribution, or `2` for a two-tailed prior distribution.
#'
#' @section Note: Please enter a positive obtained value if the observed effect
#'   is in the same directions as the effect predicted by H1, and a negative obtained value
#'   if not. This is critical to get a correct Bayes factor when you have a directional
#'   hypothesis (`tail = 1`). Note that the `sd_of_theory` must be a positive value.
#' 
#' @return The function returns the Bayes factor (evidence for H1 over H0) as a single
#' numeric value.
#' 
#' @examples 
#' # One-tailed test with the observed effect in the same direction as H1 predicted
#' SampleSizeBf::Bf_cauchy(sd = 0.4, obtained = 1, dfdata = 98, mean_of_theory = 0,
#'                         sd_of_theory = 1.5, tail = 1)
#'                         
#' # One-tailed test with the observed effect in the opposite direction as H1 predicted
#' SampleSizeBf::Bf_cauchy(sd = 0.4, obtained = -0.3, dfdata = 98, mean_of_theory = 0,
#'                         sd_of_theory = 1.5, tail = 1)
#'
#' @export
Bf_cauchy <- function(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail = c(1, 2)) {
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