calculate_df <- function(n, design = c("within", "between")) {
  # Validating arguments ---------------------------
  # n must be something
  
  if (design == "between") {
    2 * n - 2
  } else {
    n - 1
  }
}