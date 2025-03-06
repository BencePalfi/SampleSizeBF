#' Calculate sample size for Bayes factor
#' 
#' This function calculates the minimum sample size of your design
#' you may need to find evidence for H1 (or for H0). The function assumes equal 
#' group sizes for between-subject designs.
#' 
#' 
Bf_samplesize <- function(distribution_likelihood = c("normal", "t"), distribution_prior = c("normal", "cauchy", "uniform"), sd_of_theory, tail = c(1, 2), n_start, n_end, step, sd1, sd2 = NULL, threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95)) {
  
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
  
  # Get range
  if (tail == 1) {
    range <- c(0, Inf)
  } else {
    range <- c(-Inf, Inf)
  }
  
  # Get uniform min and max
  if (distribution_prior == "uniform") {
    min <- sd_of_theory - sd_of_theory
    max <- sd_of_theory + sd_of_theory
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
    ### Create likelihood dynamically
    if (distribution_likelihood == "normal") {
      data_mod <- bayesplay::likelihood(family = "normal", mean = sd_of_theory, sd = se_temp)
    } else if (distribution_likelihood == "t") {
      data_mod <- bayesplay::likelihood(family = "t", mean = sd_of_theory, sd = se_temp, df = dfdata)
    } else {
      stop("Unsupported distribution. Choose 'normal', 't'.")
    }
    
    ### Create null prior
    h0_mod <- bayesplay::prior(family = "point", point = 0)
    
    ### Create alternative prior dynamically
    if (distribution_prior == "normal") {
      h1_mod <- bayesplay::prior(family = "normal", mean = 0, sd = sd_of_theory)
    } else if (distribution_prior == "uniform") {
      h1_mod <- bayesplay::prior(family = "uniform", min = min, max = max)
    } else if (distribution_prior == "cauchy") {
      h1_mod <- bayesplay::prior(family = "cauchy", mean = 0, scale = sd_of_theory)
    } else {
      stop("Unsupported distribution. Choose 'normal', 'uniform', or 'cauchy'.")
    }
    
    m1 <- bayesplay::integral(data_mod * h1_mod)
    m0 <- bayesplay::integral(data_mod * h0_mod)
    
    bf <- m1 / m0

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
    ### Create likelihood dynamically
    if (distribution_likelihood == "normal") {
      data_mod <- bayesplay::likelihood(family = "normal", mean = 0, sd = se_temp)
    } else if (distribution_likelihood == "t") {
      data_mod <- bayesplay::likelihood(family = "t", mean = 0, sd = se_temp, df = dfdata)
    } else {
      stop("Unsupported distribution. Choose 'normal', 't'.")
    }
    
    ### Create null prior
    h0_mod <- bayesplay::prior(family = "point", point = 0)
    
    ### Create alternative prior dynamically
    if (distribution_prior == "normal") {
      h1_mod <- bayesplay::prior(family = "normal", mean = 0, sd = sd_of_theory)
    } else if (distribution_prior == "uniform") {
      h1_mod <- bayesplay::prior(family = "uniform", min = min, max = max)
    } else if (distribution_prior == "cauchy") {
      h1_mod <- bayesplay::prior(family = "cauchy", mean = 0, scale = sd_of_theory)
    } else {
      stop("Unsupported distribution. Choose 'normal', 'uniform', or 'cauchy'.")
    }

    m1 <- bayesplay::integral(data_mod * h1_mod)
    m0 <- bayesplay::integral(data_mod * h0_mod)
    
    bf <- m1 / m0
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