Bf_samplesize <- function(expected_mean, sd_of_theory, sd1, sd2 = NULL, tail = c(1, 2), threshold = c(3, 6, 10), tpr = c(.5, .8, .9, .95), design = c("between", "within"), n_start = 5, n_end = 100, step = 1) {
  ## Validating arguments
  if (design == "within" & !is.null(sd2)) cat("For within subject design sd1 will be used and sd2 will be ignored.\n")
  .tpr <- tpr
  .threshold <- threshold
  
  ## Number of subjects
  n <- seq(n_start, n_end, step)
  
  # Get stop value
  stop_value <- 
    tpr_table %>% 
    dplyr::filter(threshold == .threshold, tpr == .tpr) %>% 
    dplyr::pull(bf)

  ## Initial parameter values
  i <-  0
  bf <- 1
  n_out <- 0

  while (!(bf > stop_value | bf < 1 / stop_value)) {
    i <- i + 1

    ## Calculating standard error of the estimate based on the obtained SD
    if(design == "between") {
      SE_temp <- sqrt(((n[i] - 1) * sd1 ^ 2 + (n[i] - 1) * sd2 ^ 2) / (2 * n[i] - 2)) * sqrt(2 / n[i])
      dfdata <- 2 * n[i]-2
    } else {
      SE_temp <- sd1 / sqrt(n[i])
      dfdata <- n[i] - 1
    }

    ## Calculating the Bf with given n
    bf <- Bf(SE_temp, obtained = expected_mean, dfdata = dfdata,  meanoftheory = 0, sdtheory = sd_of_theory, tail = tail)
    n_out <- n[i]

    # Break if maximum n reached without conclusive Bf
    if (n_out == n_end) {
      stop(paste("The maximum expected Bf that you can reach with this sample size:", round(bf, 2)))
      }
    }

  cat("The expected sample size with threshold of", threshold, "and true positive rate of", tpr, "is", n_out, sep = " ")
  return(n_out)
}
