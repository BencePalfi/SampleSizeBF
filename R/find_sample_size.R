#' Determine the Required Sample Size for a Given Hypothesis Mean
#'
#' This function calculates the minimum sample size needed to obtain 
#' a conclusive Bayes factor (BF) in favor of either H1 or H0. The calculation 
#' is based on an iterative approach that increases the sample size until the 
#' Bayes factor exceeds a predefined threshold.
#'
#' @param predicted_value The mean of the hypothesis being tested. Typically, 
#'   this is either the effect size under H1 or 0 under H0.
#' @param n A numeric vector representing the sequence of sample sizes to be tested.
#' @param design A character string, either `"within"` or `"between"`, 
#'   specifying the study design.
#' @param distribution_likelihood A character string specifying the 
#'   likelihood distribution (e.g., `"normal"` or `"t"`).
#' @param distribution_prior A character string specifying the 
#'   prior distribution (e.g., `"normal"`, `"cauchy"`, or `"uniform"`).
#' @param sd1 Standard deviation of the first sample.
#' @param sd2 Standard deviation of the second sample (if applicable; `NULL` for within-subject designs).
#' @param sd_of_theory Standard deviation of the theoretical prior distribution.
#' @param stop_value The Bayes factor threshold at which the function stops iterating.
#' @param n_end The maximum sample size allowed before stopping the search.
#' @param create_likelihood A function that constructs the likelihood given 
#'   a mean, standard error, and degrees of freedom.
#' @param create_prior A function that constructs the prior distribution.
#' @param tail The number of tails (1 or 2).
#'
#' @return The smallest sample size `n` required to reach the stopping 
#'   Bayes factor threshold. Returns `NA` if the maximum sample size is 
#'   reached without reaching conclusive evidence.
#'
#' @details
#' The function iterates through increasing sample sizes, calculating the 
#' corresponding standard error and degrees of freedom. For each sample size, 
#' it computes the Bayes factor using the provided likelihood and prior 
#' distributions. The process stops when the Bayes factor exceeds the 
#' predefined threshold (`stop_value`) or when the maximum allowed sample 
#' size (`n_end`) is reached.
#'
#' @seealso \code{\link{Bf_samplesize}}, \code{\link{bayesplay::likelihood}}, 
#'   \code{\link{bayesplay::prior}}
#'
#' @examples
#' \dontrun{
#' # Example: Find required sample size for a given hypothesis mean
#' required_n <- find_sample_size(
#'   predicted_value = 0.5,
#'   n = seq(5, 100, 5),
#'   design = "within",
#'   distribution_likelihood = "normal",
#'   distribution_prior = "cauchy",
#'   sd1 = 1,
#'   sd2 = NULL,
#'   sd_of_theory = 0.5,
#'   stop_value = 6,
#'   n_end = 100,
#'   create_likelihood = create_likelihood,
#'   create_prior = create_prior
#' )
#' print(required_n)
#' }
find_sample_size <- function(
    predicted_value, 
    n, 
    design, 
    distribution_likelihood, 
    distribution_prior, 
    sd1, 
    sd2, 
    sd_of_theory, 
    stop_value, 
    n_end, 
    create_likelihood, 
    create_prior,
    tail
) {
  i <- 0
  bf <- 1
  n_out <- NA_integer_

  while (!(bf > stop_value | bf < 1 / stop_value)) {
    i <- i + 1
    
    if (i > length(n)) break  # Stop if we exceed the sample size range
    
    # Calculate standard error
    se_temp <- if (design == "between") {
      calculate_se(sd1 = sd1, sd2 = sd2, n = n[i])
    } else {
      calculate_se(sd1 = sd1, n = n[i])
    }
    
    dfdata <- if (design == "between") 2 * n[i] - 2 else n[i] - 1
    
    # Generate likelihood and priors
    data_mod <- create_likelihood(distribution_likelihood, mean = predicted_value, se = se_temp, dfdata = dfdata)
    h0_mod <- bayesplay::prior(family = "point", point = 0)
    h1_mod <- create_prior(distribution_prior, sd_of_theory, tail = tail)
    
    # Compute Bayes factor
    m1 <- bayesplay::integral(data_mod * h1_mod)
    m0 <- bayesplay::integral(data_mod * h0_mod)
    bf <- m1 / m0

    n_out <- n[i]
    
    # Stop if maximum sample size is reached without conclusive BF
    if (n_out == n_end) {
      n_out <- NA_integer_
      break
    }
  }
  
  return(n_out)
}
