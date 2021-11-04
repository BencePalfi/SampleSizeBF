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

#' title
#' 
#' The function assumes equal group sizes for between-subject designs.
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
      h1_bf = h1_bf,
      h0 = h0,
      h0_bf = h0_bf)
    )
}
