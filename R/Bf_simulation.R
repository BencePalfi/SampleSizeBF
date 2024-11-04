#' Function simulating studies and Bayesian sequential analyses
#' 
#' This function...
#' 
#' @param sd_of_theory numeric.
#' @param sd1 numeric.
#' @param sd2 numeric.
#' @param correlation numeric.
#' @param threshold integer.
#' @param stopping_rule character.
#' @param n integer.
#' @param iterations integer.
#' @param Bf_calculation function object.
#' @param tail integer.
#' @param true_effect logical.
#' 
#' @return The function returns a tibble with two columns: `id` that contains
#'   the unique identifier of the given iteration, and `hit` which contains whether
#'   the resulting Bayes factor was larger than the specified threshold (`2`),
#'   inconclusive (`0`), or smaller than 1 / threshold (`1`).
#' 
#' @export
simulate_Bf <- function(
  sd_of_theory,
  sd1,
  sd2 = NULL,
  correlation = 0,
  threshold = c(3, 6, 10),
  stopping_rule = c("optional", "fixed"),
  n = 100,
  iterations = 5,
  Bf_calculation = c("Bf_normal", "Bf_cauchy"),
  tail = c(1, 2),
  true_effect = TRUE) {
  # Choose the appropriate function based on Bf_type
  Bf_calculation_func <- switch(Bf_calculation,
                           "normal" = Bf_normal,
                           "cauchy" = Bf_cauchy)
  
  # Progress bar
  pb <- progress::progress_bar$new(
      format = " simulation progress [:bar] :percent eta: :eta",
      total = iterations, clear = FALSE, width = 60)
  
  # Calculate results
  res <- 
    tibble::tibble(
      id = 1:iterations,
      hit = purrr::map_int(1:iterations, 
                       function(x) {
                         pb$tick() # progress bar tick (useful to track progress of long analyses)
                         
                         sample <- generate_sample(
                           n = n,
                           sd_of_theory = sd_of_theory,
                           sd1 = sd1, 
                           sd2 = sd2,
                           correlation = correlation,
                           true_effect = true_effect
                         )
                         
                         if (stopping_rule == "fixed") {
                           bf <- fixed_BF(
                             sample = sample, 
                             sd_of_theory = sd_of_theory,
                             tail = tail,
                             Bf_calculation = Bf_calculation_func)
                         } else if (stopping_rule == "optional") {
                           bf <- optional_BF(
                             sample = sample,
                             threshold = threshold,
                             sd_of_theory = sd_of_theory,
                             tail = tail,
                             Bf_calculation = Bf_calculation_func)
                         }
                         
                         if (bf >= threshold) {
                           return(2L)
                         } else if (bf <= 1 / threshold) {
                           return(1L)
                         } else {
                           return(0L)
                         }
                       })
    )
  
  return(res)
}

#' Function for generating samples
#' 
#' This function generates a sample with a given sample size for a true effect (with a given hypothetical effect size) or a true null effect.
#' The function can generate samples for independent and paired study designs. For independent study designs provide both sd1 (experimental group) and 
#' s2 (control group), while for paired study designs provide only sd1 for the standard deviation of the difference scores between the two groups. For
#' independent study designs the sample generating function assumes 0 correlation between the two groups.
#' 
#' @param n integer.
#' @param sd_of_theory numeric.
#' @param sd1 numeric.
#' @param sd2 numeric.
#' @param true_effect logical.
#' @param correlation numeric.
#' 
#' @return The function returns a tibble with ...
#' 
#' @export
generate_sample <- function(n, sd_of_theory = 0, sd1, sd2 = NULL, true_effect = TRUE, correlation = 0) {
  if (true_effect) {
    true_effect_size <- sd_of_theory
  } else {
    true_effect_size <- 0
  }
  
  if (is.null(sd2)) {
    sample <- tibble::tibble(
      outcome_diff = rnorm(mean = 0, sd = 1, n = n)
    ) %>% 
      dplyr::mutate(
        outcome_diff = outcome_diff * sd1 + true_effect_size 
      )
  } else {
    # generate two correlated random variables with mean = 0 and sd = 1
    sample <- as.data.frame(
      MASS::mvrnorm(
        n,
        mu = c(0, 0),
        Sigma = matrix(
            c(1, correlation, correlation, 1),
          ncol = 2,
          byrow = TRUE),
        empirical = F)) %>% 
      dplyr::rename(
        outcome_exp = V1,
        outcome_con = V2
      ) %>% 
      # Transform these variables to match mean and sd from pilot study / assumed parameters
      dplyr::mutate(
        outcome_exp = outcome_exp * sd1 + true_effect_size,
        outcome_con = outcome_con * sd2
      )
  }
  
  return(sample)
}

#' Function to calculate Bayes factor with fixed stopping
#' 
#' This function ...
#' 
#' @param sample tibble.
#' @param sd_of_theory numeric.
#' @param tail integer.
#' @param Bf_calculation function object.
#' 
#' @return The function returns
#' 
#' @export
fixed_BF <- function (sample, sd_of_theory, tail = c(1, 2), Bf_calculation) {
  # Calculate the t test for simulated data (one sample only)
  if (length(sample) != 1) {
    t_test <- t.test(sample$outcome_exp, sample$outcome_con, paired = FALSE, var.equal = TRUE) %>% broom::tidy()
    } else {
      t_test <- t.test(sample$outcome_diff, mu = 0) %>% broom::tidy()
    }
  
  # Calculate Bayes factor and other parameters
  bf <- Bf_calculation(
    sd = abs(t_test$estimate / t_test$statistic),
    obtained = t_test$estimate,
    dfdata = t_test$parameter,
    mean_of_theory = 0,
    sd_of_theory = sd_of_theory,
    tail = tail
    )
  
  attributes(bf) <- NULL
  
  return(bf)
}

#' Function to calculate Bayes factor with optional stopping
#' 
#' The function ...
#' 
#' @param sample tibble.
#' @param threshold integer.
#' @param sd_of_theory numeric.
#' @param tail integer.
#' @param Bf_calculation function object.
#' 
#' @return The function returns
#' 
#' @export
optional_BF <- function(sample, threshold, sd_of_theory, tail = c(1, 2), Bf_calculation) {
  i <-  0
  bf <- 1
  n_out <- 0
  n <- 5 : nrow(sample)
  
  while (!(bf > threshold | bf < 1 / threshold)) { 
    i <- i + 1
    
    # Get a slice of the sample
    data_current <- dplyr::slice(sample, 1 : n[i])
    
    # Calculate the t test for simulated data (one sample only)
    if (length(data_current) != 1) {
      t_test <- t.test(data_current$outcome_exp, data_current$outcome_con, paired = FALSE, var.equal = TRUE) %>% broom::tidy()
      } else {
        t_test <- t.test(data_current$outcome_diff, mu = 0) %>% broom::tidy()
      }
    
    # Calculate Bayes factor and other parameters
    bf <- Bf_calculation(
      sd = abs(t_test$estimate / t_test$statistic),
      obtained = t_test$estimate,
      dfdata = t_test$parameter, 
      mean_of_theory = 0,
      sd_of_theory = sd_of_theory,
      tail = tail)
    
    attributes(bf) <- NULL
    
    n_out <- n[i]
    
    # Break if maximum n reached without conclusive Bf
    if (n_out == max(n)) {
      return(bf)
    }
  }
  
  return(bf)
}

