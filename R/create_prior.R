#' Generate a Prior Distribution for Bayes Factor Calculation
#'
#' This function creates a prior distribution for use in Bayesian inference 
#' based on the specified distribution family.
#'
#' @param distribution_prior A character string specifying the prior distribution. 
#'   Options include `"normal"`, `"uniform"`, and `"cauchy"`.
#' @param sd_of_theory The standard deviation or scale parameter for the prior 
#'   (depends on the chosen distribution).
#' @param tail The number of tails (1 or 2).
#'
#' @return A `bayesplay` prior object corresponding to the specified distribution.
#'
#' @details
#' - **Normal prior:** Centered at `0` with standard deviation `sd_of_theory`.
#' - **Uniform prior:** A flat distribution ranging from `-sd_of_theory` to `+sd_of_theory`.
#' - **Cauchy prior:** Centered at `0` with scale `sd_of_theory`.
#'
#' @seealso \code{\link{bayesplay::prior}}, \code{\link{Bf_samplesize}}
#'
#' @examples
#' \dontrun{
#' # Normal prior with SD = 0.5
#' normal_prior <- create_prior("normal", sd_of_theory = 0.5)
#' 
#' # Uniform prior from -1 to 1
#' uniform_prior <- create_prior("uniform", sd_of_theory = 1)
#'
#' # Cauchy prior with scale = 0.3
#' cauchy_prior <- create_prior("cauchy", sd_of_theory = 0.3)
#' }
create_prior <- function(distribution_prior, sd_of_theory, tail) {
  if (distribution_prior == "normal") {
    range <- if (tail == 1) c(0, Inf) else c(-Inf, Inf)
    return(bayesplay::prior(family = "normal", mean = 0, sd = sd_of_theory, range = range))
  } else if (distribution_prior == "uniform") {
    if (tail == 1) {
      min <- 0
      max <- 2 * sd_of_theory
    } else {
      min <- -2 * sd_of_theory
      max <-  2 * sd_of_theory
    }
    return(bayesplay::prior(family = "uniform", min = min, max = max))
  } else if (distribution_prior == "cauchy") {
    range <- if (tail == 1) c(0, Inf) else c(-Inf, Inf)
    return(bayesplay::prior(family = "cauchy", location = 0, scale = sd_of_theory, range = range))
  } else {
    stop("Unsupported prior distribution. Choose 'normal', 'uniform', or 'cauchy'.")
  }
}
