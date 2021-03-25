#' Function simulating studies and Bayesian sequential analyses
#' 
#' This function...
#' 
#' @param N sample size
#' @param mean_exp Mean in the experimental condition
#' 
#' @return The function returns four values
#' 
#' @export
simulate_BF <- function(N,
                        mean_exp,
                        SD_exp,
                        mean_con,
                        SD_con,
                        correlation = 0,
                        sdtheory,
                        cutoff = 3,
                        stopping_rule = c('optional','fixed'),
                        steps = 10,
                        iterations = 1000){
  
  pb <- progress::progress_bar$new(
      format = " simulation progress [:bar] :percent eta: :eta",
      total = iterations, clear = FALSE, width= 60)
  
  pb$tick() # progress bar tick (useful to track progress of long analyses)
  
  # generate two correlated random variables with mean = 0 and sd = 1
  data = as.data.frame(MASS::mvrnorm(N,
                                     mu = c(0, 0),
                                     Sigma = matrix(
                                       c(1, correlation, correlation, 1),
                                       ncol = 2,
                                       byrow = TRUE
                                       ),
                                     empirical = F))
  names(data) = c("outcome_exp", "outcome_con")
  
  # transform these variables to match mean and sd from pilot study / assumed parameters
  data[,"outcome_exp"] = data[,"outcome_exp"] * SD_exp + mean_exp
  data[,"outcome_con"] = data[,"outcome_con"] * SD_con + mean_con
  
  #Calculating the BFs
  if (stopping_rule == 'optional') {
    for(i in seq(10,N,steps)){ 
      data_current = data[1:i,]
      effect <- t.test(data_current[,"outcome_exp"], data_current[,"outcome_con"], paired=FALSE, var.equal = TRUE)
      mean_dif = as.numeric(effect$estimate[1]) - as.numeric(effect$estimate[2])
      se_dif = abs(as.numeric(mean_dif/effect$statistic[1]))
      BF = Bf(sd = se_dif, obtained = mean_dif, dfdata = effect$parameter,  meanoftheory = 0, sdtheory = sdtheory, tail = 1)
      if(BF > cutoff){break}
      if(BF < 1/cutoff){break}
    }
    sample_size = nrow(data_current)
  }
  if (stopping_rule == 'fixed') {
    effect <- t.test(data[,"outcome_exp"], data[,"outcome_con"], paired=FALSE, var.equal = TRUE)
    mean_dif = as.numeric(effect$estimate[1]) - as.numeric(effect$estimate[2])
    se_dif = abs(as.numeric(mean_dif/effect$statistic[1]))
    BF = Bf(sd = se_dif, obtained = mean_dif, dfdata = effect$parameter,  meanoftheory = 0, sdtheory = sdtheory, tail = 1)
    t_value = as.numeric(effect$statistic[1])
    sample_size = nrow(data)
  }
  
  cor = cor(data[,"outcome_exp"], data[,"outcome_con"])
  
  # the function returns the Bayes Factor value at the stopping point of the study
  return(c(BF, mean_dif, se_dif,cor, sample_size))
}