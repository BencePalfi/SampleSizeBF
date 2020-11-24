Bf_samplesize <- function(n_start = 5, n_end = 100, step = 1, design = c("between", "within"), sd1, sd2 = NULL, expected_effect, tail, threshold = 3) {
  ##number of subjects
  n <- seq(n_start, n_end, step)
  i = 0
  x = 1
  n_out = 0
  while(x < threshold| x > 1/threshold){
    i = i + 1
    ##standard error of the estimate based on the obtained SD
    if(design == "between"){
      SE_temp <- sqrt( ((n[i]-1)*sd1^2 + (n[i]-1)*sd2^2) /(2*n[i]-2))*sqrt(2/n[i])
      dfdata = 2*n[i]-2
    } else {
      SE_temp <- sd1/sqrt(n[i])
      dfdata= n[i]-1
    } 
    
    ##B
    x <- Bf(SE_temp, obtained = expected_effect, dfdata = dfdata,  meanoftheory = 0, sdtheory = expected_effect, tail= tail)
    n_out <- n[i]
    if(n[i]==n_end) {
      cat("The maximum expected Bf that you can reach with this sample size:", x[i], sep = " ")
      return(x[i])
      break
    }
  }
  cat("The expected sample size with", x, "is", n_out,sep = " ")
  return(n_out)
  }