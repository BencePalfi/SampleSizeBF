calculate_se <- function(sd1, sd2 = NULL, n) {
  # Validating arguments ---------------------------
  #sd1 must be something
  
  if (is.null(sd2)) {
    sd1 / sqrt(n)
  } else {
    sqrt(((n - 1) * sd1 ^ 2 + (n - 1) * sd2 ^ 2) / (2 * n - 2)) * sqrt(2 / n)
    }
}