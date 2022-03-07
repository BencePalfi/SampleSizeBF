#' Calculate sample size for Bayes factor
#' 
#' This function calculates the minimum sample size of your design
#' you may need to find evidence for H1 (or for H0). The function assumes equal 
#' group sizes for between-subject designs.
#' 
#' 
Bf_samplesize <- function(Bf_calculation, n_start, n_end, step, sd1, sd2 = NULL, threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95), ...) {
  call_parameters <- rlang::dots_list(...)
  
  # Validating arguments ---------------------------
  if (n_start < 5) stop("The minimum sample size should be at least 5.")
  .tpr <- tpr
  .threshold <- threshold
  # sd_of_theory should not be null
  # sd1 must be something
  # step must be integer and minimum 1
  # tail must be 1 or 2
  # threshold must be integer and either 3, 6, 10
  
  # Number of subjects ---------------------------
  n <- seq(n_start, n_end, step)
  
  # Get stop value ---------------------------
  stop_value <-
    tpr_table %>%
    dplyr::filter(threshold == .threshold, tpr == .tpr) %>%
    dplyr::pull(bf)
  
  # Get design ---------------------------
  if (is.null(sd2)) {
    design <- "within"
  } else {
    design <- "between"
  }
  
  # Initial parameter values ---------------------------
  i <-  0
  bf <- 1
  n_out <- 0
  
  # Sample size for H1 ---------------------------
  while (!(bf > stop_value | bf < 1 / stop_value)) {
    i <- i + 1
    
    ## Calculating standard error of the estimate based on the obtained SD
    # Between design
    if(design == "between") {
      se_temp <- calculate_se(sd1 = sd1, sd2 = sd2, n = n[i])
      dfdata <- 2 * n[i] - 2
    } else {
      se_temp <- calculate_se(sd1 = sd1, n = n[i])
      dfdata <- n[i] - 1
    }
    
    ## Calculating the Bf with given n
    bf_quote <- rlang::enexpr(Bf_calculation)
    bf_data <- rlang::dots_list(sd = se_temp, dfdata = dfdata, !!!call_parameters)
    bf <- rlang::eval_tidy(bf_quote, data = bf_data)
    # bf <- Bf_calculation(sd = se_temp, dfdata = dfdata, ...)
    n_out <- n[i]
    
    # Break if maximum n reached without conclusive Bf
    if (n_out == n_end) {
      n_out <- NA_integer_
      break
    }
  }
  
  h1 <- as.integer(n_out)
  h1_bf <- bf
  
  
  # Reset parameter values ---------------------------
  i <-  0
  bf <- 1
  n_out <- 0
  obtained <- 0
  
  # Sample size for H0 ---------------------------
  while (!(bf > stop_value | bf < 1 / stop_value)) {
    i <- i + 1
    
    ## Calculating standard error of the estimate based on the obtained SD
    if(design == "between") {
      se_temp <- calculate_se(sd1 = sd1, sd2 = sd2, n = n[i])
      dfdata <- 2 * n[i] - 2
    } else {
      se_temp <- calculate_se(sd1 = sd1, n = n[i])
      dfdata <- n[i] - 1
    }
    
    ## Calculating the Bf with given n
    bf_data <- rlang::dots_list(sd = se_temp, dfdata = dfdata, obtained = 0, !!!call_parameters, .homonyms = "first")
    bf <- rlang::eval_tidy(bf_quote, data = bf_data)
    n_out <- n[i]
    
    # Break if maximum n reached without conclusive Bf
    if (n_out == n_end) {
      n_out <- NA_integer_
      break
    }
  }
  
  h0 <- as.integer(n_out)
  h0_bf <- bf
  
  # Return output ---------------------------
  if (design == "between") {
    cat("The expected sample size per group with threshold of", threshold, "and true positive rate of", tpr, "is", "H1:", h1, "H0:", h0, "\n", sep = " ")
  } else {
    cat("The expected sample size with threshold of", threshold, "and true positive rate of", tpr, "is", "H1:", h1, "H0:", h0, "\n", sep = " ")
  }
  
  if (is.na(h1) | is.na(h0)) {message("Note: Missing values indicate that the maximum available sample size is not sufficient to reach the desired threshold.")}
  
  invisible(
    list(
      h1 = h1,
      h0 = h0)
  )
}

#'  Calculate the sample size for Bayes factor with normal prior
#'  
#' This function calculates the minimum sample size of your design
#' you may need to find evidence for H1 (or for H0). The function
#' assumes that the model of H1 is a good representation of the real
#' effect. The function models the predictions of H1 with a normal
#' distribution, and assumes equal group sizes for between-subject designs.
#'  
#'  @param sd_of_theory numeric. The standard deviation of the model of H1 (i.e. the expected effect size).
#'  @param sd1 numeric. The expected standard deviation of group1, or the standard deviation of the difference scores in case of a whithin-subject design.
#'  @param sd2 numeric. The expected standard deviation of group2, or `NULL` in case of a within-subject design.
#'  @param tail integer. Either `1` for a half-normal prior distribution, or `2` for a two-tailed prior distribution.
#'  @param threshold integer. The required level of evidence for H1. Either `3`, `6`, or `10`.
#'  @param tpr numeric. The required true positive rate. The long-run probability of finding a sensitive Bayes factor (cf. power in NHST). Either `.5`, `.8`, `.9`, or `.95`.
#'  @param n_start integer. The minimum sample size of interest for one group. Must be at least `5`.
#'  @param n_end integer. The maximum sample size of interest for one group.
#'  @param step integer. Divides the range of sample sizes of interest between `n_start`and `n_end`. The smallest values is 1.
#'  
#'  @section Note: The runtime of the function can increase substantially by increasing
#'    the `threshold`, `tpr`, and `n_end`. If the runtime takes too long consider increasing
#'    the `step` parameter. Bare in mind that increasing `step` may produce an inflated
#'    minimum sample size, thus it is up to you to find the balance between precision and
#'    computational speed.
#'  
#'  @return The function returns a list of 2 numeric vectors. `h1` is the minimum required
#'    sample size that your design needs to reach a sensitive Bayes factor with the given
#'    true positive rate assuming that H1 is true. `h0` is the minimum required
#'    sample size that your design needs to reach a sensitive Bayes factor with the given
#'    true positive rate assuming that H0 is true.
#' 
#' @export
Bf_samplesize_normal <- function(sd_of_theory, sd1, sd2 = NULL, tail = c(1, 2), threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95), n_start = 5, n_end = 100, step = 1) {
  Bf_samplesize(Bf_normal(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail),
                sd1 = sd1,
                sd2 = sd2,
                mean_of_theory = 0,
                sd_of_theory = sd_of_theory,
                obtained = sd_of_theory,
                tail = tail,
                threshold = threshold,
                tpr = tpr,
                n_start = n_start,
                n_end = n_end,
                step = step
                )
}

#'  Function to calculate the samplesize for Bayes factor with uniform prior
#'  
#'  This function ...
#'  
#'  @param upper numeric.
#'  @param sd1 numeric.
#'  @param sd2 numeric.
#'  @param tail integer.
#'  @param threshold integer.
#'  @param tpr numeric.
#'  @param n_start integer.
#'  @param n_end integer.
#'  @param step integer.
#'  
#'  @return The function returns ...
#' 
#' @export
Bf_samplesize_uniform <- function(upper, sd1, sd2 = NULL, tail = c(1, 2), threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95), n_start = 5, n_end = 100, step = 1) {
  # Validation
  # upper must be positive
  
  if (tail == 1) {
    lower <- 0
  } else {
    lower <- -upper
  }
  
  Bf_samplesize(Bf_uniform(sd, obtained, dfdata, lower, upper),
                sd1 = sd1,
                sd2 = sd2,
                upper = upper,
                lower = lower,
                obtained = upper / 2,
                tail = tail,
                threshold = threshold,
                tpr = tpr,
                n_start = n_start,
                n_end = n_end,
                step = step
                )
}

#'  Function to calculate the samplesize with Cauchy prior
#'  
#'  This function ...
#'  
#'  @param sd_of_theory numeric.
#'  @param sd1 numeric.
#'  @param sd2 numeric.
#'  @param tail integer.
#'  @param threshold integer.
#'  @param tpr numeric.
#'  @param n_start integer.
#'  @param n_end integer.
#'  @param step integer.
#'  
#'  @return The function returns ...
#' 
#' @export
Bf_samplesize_cauchy <- function(sd_of_theory, sd1, sd2 = NULL, tail = c(1, 2), threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95), n_start = 5, n_end = 100, step = 1) {
  Bf_samplesize(Bf_cauchy(sd, obtained, dfdata, mean_of_theory, sd_of_theory, tail),
                sd1 = sd1,
                sd2 = sd2,
                mean_of_theory = 0,
                obtained = sd_of_theory,
                sd_of_theory = sd_of_theory,
                tail = tail,
                threshold = threshold,
                tpr = tpr,
                n_start = n_start,
                n_end = n_end,
                step = step
                )
}