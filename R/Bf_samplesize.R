Bf_samplesize <- function(n_start = 5, n_end = 100, step = 1, design = c("between", "within"), sd1, sd2 = NULL, expected_effect, tail = 1, threshold = 3) {
  ## Validating arguments
  if (design == "within" & !is.null(sd2)) cat("For within subject design sd1 will be used and sd2 will be ignored.\n")
  
  ## Number of subjects
  n <- seq(n_start, n_end, step)
  
  ## Initial parameter values
  i <-  0
  bf <- 1
  n_out <- 0
  
  while (!(bf > threshold | bf < 1 / threshold)) {
    i <- i + 1
    
    ## Calculating standard error of the estimate based on the obtained SD
    if(design == "between"){
      SE_temp <- sqrt(((n[i] - 1) * sd1 ^ 2 + (n[i] - 1) * sd2 ^ 2) / (2 * n[i] - 2)) * sqrt(2 / n[i])
      dfdata <- 2 * n[i]-2
    } else {
      SE_temp <- sd1 / sqrt(n[i])
      dfdata <- n[i] - 1
    } 
    
    ## Calculating the Bf with given n
    bf <- Bf(SE_temp, obtained = expected_effect, dfdata = dfdata,  meanoftheory = 0, sdtheory = expected_effect, tail = tail)
    n_out <- n[i]
    
    # Break if maximum n reached without conclusive Bf
    if (n[i] == n_end) {
      cat("The maximum expected Bf that you can reach with this sample size:", bf, sep = " ")
      return(bf)
      break
      }
    }
  
  cat("The expected sample size with", bf, "is", n_out, sep = " ")
  return(n_out)
}
