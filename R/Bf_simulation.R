#' Function simulating studies and Bayesian sequential analyses
#' 
#' This function...
#' 
#' @param mean_exp Mean in the experimental condition
#' 
#' @return The function returns four values
#' 
#' @export
simulate_Bf <- function(
  sd_of_theory,
  sd1,
  sd2 = NULL,
  correlation = 0,
  threshold = 3,
  stopping_rule = c("optional", "fixed"),
  n = 100,
  iterations = 5,
  Bf_calculation = c(Bf_normal, Bf_cauchy),
  tail = c(1, 2)) {
  
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
                           correlation = correlation
                         )
                         
                         if (stopping_rule == "fixed") {
                           bf <- fixed_BF(
                             sample = sample, 
                             sd_of_theory = sd_of_theory,
                             tail = tail,
                             Bf_calculation = Bf_calculation)
                         } else if (stopping_rule == "optional") {
                           bf <- optional_BF(
                             sample = sample,
                             threshold = threshold,
                             sd_of_theory = sd_of_theory,
                             tail = tail,
                             Bf_calculation = Bf_calculation)
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

generate_sample <- function(n, sd_of_theory, sd1, sd2 = NULL, correlation = correlation) {
  if (is.null(sd2)) {
    sample <- tibble::tibble(
      outcome_diff = rnorm(mean = 0, sd = 1, n = n)
    ) %>% 
      dplyr::mutate(
        outcome_diff = outcome_diff * sd1 + sd_of_theory 
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
        outcome_exp = outcome_exp * sd1 + sd_of_theory,
        outcome_con = outcome_con * sd2
      )
  }
  
  return(sample)
}

fixed_BF <- function (sample, sd_of_theory, tail = c(1, 2), Bf_calculation) {
  # Calculate the t test for simulated data (one sample only)
  if (length(sample) != 1) {
    t_test <- t.test(sample$outcome_exp, sample$outcome_con, paired = FALSE) %>% broom::tidy()
    } else {
      t_test <- t.test(sample$outcome_diff, mu = 0, paired = FALSE, var.equal = TRUE) %>% broom::tidy()
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
        t_test <- t.test(data_current$outcome_diff, mu = 0, paired = FALSE, var.equal = TRUE) %>% broom::tidy()
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

