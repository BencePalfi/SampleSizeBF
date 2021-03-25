## TPR - threshold table
library(tidyverse)

tpr_table <- 
  tibble::tibble(
    threshold = rep(c(3, 6, 10), 4),
    tpr = c(rep(.5, 3), rep(.8, 3), rep(.9, 3), rep(.95, 3)),
    bf = c(3, 6, 10, 20, 40, 85, 70, 150, 350, 220, 520, 1370)
  )
  
usethis::use_data(tpr_table, overwrite = TRUE)
