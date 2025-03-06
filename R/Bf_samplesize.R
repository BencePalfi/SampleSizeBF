#' Calculate Sample Size for Bayes Factor
#'
#' This function calculates the minimum sample size required to obtain evidence 
#' for H1 or H0 using Bayes factors. It supports different likelihood and prior 
#' distributions and assumes equal group sizes for between-subject designs.
#' 
#' @param distribution_likelihood Distribution for the likelihood function ("normal" or "t").
#' @param distribution_prior Distribution for the prior ("normal", "cauchy", or "uniform").
#' @param sd_of_theory The standard deviation of the theoretical distribution.
#' @param tail The number of tails (1 or 2).
#' @param n_start Starting sample size.
#' @param n_end Maximum sample size.
#' @param step Step size for sample size increments.
#' @param sd1 Standard deviation of the first sample.
#' @param sd2 Standard deviation of the second sample (for between-subject designs).
#' @param threshold Bayes factor threshold for stopping.
#' @param tpr True positive rate.
#' 
#' @return A list with required sample sizes for H1 and H0.
#' @export
Bf_samplesize <- function(
    distribution_likelihood = c("normal", "t"), 
    distribution_prior = c("normal", "cauchy", "uniform"), 
    sd_of_theory, 
    tail = c(1, 2), 
    n_start, 
    n_end, 
    step, 
    sd1, 
    sd2 = NULL, 
    threshold = c(3, 6, 10), 
    tpr = c(.5, .8, .9, .95)
) {
  
  .tpr <- tpr
  .threshold <- threshold
  
  if (n_start < 5) stop("The minimum sample size should be at least 5.")
  if (!is.numeric(sd_of_theory) || sd_of_theory <= 0) stop("sd_of_theory must be a positive number.")
  
  # Determine the study design
  design <- if (is.null(sd2)) "within" else "between"
  
  # Generate sample sizes
  n <- seq(n_start, n_end, step)
  
  # Get the stopping Bayes factor value
  stop_value <- dplyr::filter(tpr_table, threshold == .threshold, tpr == .tpr) %>%
    dplyr::pull(bf)
  
  # Compute sample sizes for H1 and H0 ---------------------------
  h1 <- find_sample_size(
    predicted_value = sd_of_theory,
    n = n,
    design = design,
    distribution_likelihood = distribution_likelihood,
    distribution_prior = distribution_prior,
    sd1 = sd1,
    sd2 = sd2,
    sd_of_theory = sd_of_theory,
    stop_value = stop_value,
    n_end = n_end,
    create_likelihood = create_likelihood,
    create_prior = create_prior,
    tail = tail
  )
  h0 <- find_sample_size(
    predicted_value = 0,
    n = n,
    design = design,
    distribution_likelihood = distribution_likelihood,
    distribution_prior = distribution_prior,
    sd1 = sd1,
    sd2 = sd2,
    sd_of_theory = sd_of_theory,
    stop_value = stop_value,
    n_end = n_end,
    create_likelihood = create_likelihood,
    create_prior = create_prior,
    tail = tail
  )
  
  # Return output ---------------------------
  message_text <- if (design == "between") {
    sprintf("The expected sample size per group with threshold of %s and true positive rate of %s is H1: %s H0: %s", 
            threshold, tpr, h1, h0)
  } else {
    sprintf("The expected sample size with threshold of %s and true positive rate of %s is H1: %s H0: %s", 
            threshold, tpr, h1, h0)
  }
  
  cat(message_text, "\n")
  
  if (is.na(h1) || is.na(h0)) {
    message("Note: Missing values indicate that the maximum available sample size is not sufficient to reach the desired threshold.")
  }
  
  invisible(list(h1 = h1, h0 = h0))
}