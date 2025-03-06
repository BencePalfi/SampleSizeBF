#' Generate a Likelihood Distribution for Bayes Factor Calculation
#'
#' This function creates a likelihood distribution based on the specified 
#' distribution family, mean, standard error, and degrees of freedom.
#'
#' @param distribution_likelihood A character string specifying the likelihood 
#'   distribution. Options include `"normal"` and `"t"`.
#' @param mean The mean of the likelihood distribution.
#' @param se The standard error of the likelihood distribution.
#' @param dfdata The degrees of freedom (only required for `"t"` distribution).
#'
#' @return A `bayesplay` likelihood object corresponding to the specified distribution.
#'
#' @details
#' - **Normal likelihood:** Uses `mean` and `sd = se`.
#' - **T-distribution likelihood:** Uses `mean`, `sd = se`, and `df = dfdata`.
#'
#' @seealso \code{\link{bayesplay::likelihood}}, \code{\link{Bf_samplesize}}
#'
#' @examples
#' \dontrun{
#' # Normal likelihood with mean = 0 and SE = 1
#' normal_likelihood <- create_likelihood("normal", mean = 0, se = 1)
#' 
#' # T-distribution likelihood with mean = 0, SE = 1, df = 10
#' t_likelihood <- create_likelihood("t", mean = 0, se = 1, dfdata = 10)
#' }
create_likelihood <- function(distribution_likelihood, mean, se, dfdata = NULL) {
  if (distribution_likelihood == "normal") {
    return(bayesplay::likelihood(family = "normal", mean = mean, sd = se))
  } else if (distribution_likelihood == "t") {
    if (is.null(dfdata)) stop("For 't' distribution, degrees of freedom (dfdata) must be provided.")
    return(bayesplay::likelihood(family = "student_t", mean = mean, sd = se, df = dfdata))
  } else {
    stop("Unsupported likelihood distribution. Choose 'normal' or 't'.")
  }
}
