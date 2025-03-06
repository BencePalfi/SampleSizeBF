# Library ----------------------------------------------------------------------
library(tidyverse)
library(future.apply)
library(here)
devtools::load_all()

# Functions --------------------------------------------------------------------
# Create a function to safely run the simulation
safe_simulate_Bf <- purrr::safely(simulate_Bf)
Bf_samplesize_normal_safe <- purrr::safely(Bf_samplesize_normal)
Bf_samplesize_cauchy_safe <- purrr::safely(Bf_samplesize_cauchy)

# Run calculations -------------------------------------------------------------
# Calculation configurations ---------------------------------------------------
# Dataframe of the input parameters for the simulation
simulation_options <-
  expand.grid(
    sd_of_theory = 30,
    sd1 = 87.34,
    tpr = c(.5, .8, .9, .95),
    threshold = c(3, 6, 10),
    stopping_rule = c("fixed", "optional"),
    Bf_calculation = c("normal", "cauchy"),
    tail = c(1, 2),
    true_effect = c(TRUE, FALSE)
  ) |> 
  as_tibble() |> 
  mutate(options_id = row_number()) |> 
  # Nest all iterations as a list
  nest(input = -options_id)

# Calculate samplesize with the input parameters
simulation_samplesize <- 
  simulation_options |> 
  mutate(samplesize = map(
    input,
    ~ {
      if (.x$Bf_calculation == "normal") {
        Bf_samplesize_normal_safe(
          sd_of_theory = 30,
          sd1 = 87.34,
          sd2 = NULL,
          tail = .x$tail,
          threshold = .x$threshold,
          tpr = .x$tpr,
          n_start = 5,
          n_end = 10000,
          step = 1
        )
      } else if (.x$Bf_calculation == "cauchy") {
        Bf_samplesize_cauchy_safe(
          sd_of_theory = 30,
          sd1 = 87.34,
          sd2 = NULL,
          tail = .x$tail,
          threshold = .x$threshold,
          tpr = .x$tpr,
          n_start = 5,
          n_end = 10000,
          step = 1
        )
      } else {
        stop("Unknown Bf_calculation value")
      }
    }
  ))

# Extract the sample sizes
simulation_samplesize <-
  simulation_samplesize |>
  dplyr::mutate(
    h1 = purrr::map_dbl(samplesize, ~ { if (!is.null(.x$error)) NA_real_ else .x$result$h1 }),
    h0 = purrr::map_dbl(samplesize, ~ { if (!is.null(.x$error)) NA_real_ else .x$result$h0 })
  ) |>
  tidyr::pivot_longer(
    cols = c(h1, h0),
    names_to = "hypothesis",
    values_to = "n"
  )

usethis::use_data(simulation_samplesize, overwrite = TRUE)

# Run simulations
simulation_outputs <-
  simulation_samplesize |>
  mutate(output = map2(
    input,
    n,
    ~ {
      if (is.na(.y)) {
        NULL
      } else {
        safe_simulate_Bf(
          sd_of_theory = .x$sd_of_theory,
          sd1 = .x$sd1,
          sd2 = NULL,
          correlation = 0,
          threshold = .x$threshold,
          stopping_rule = .x$stopping_rule,
          n = .y, 
          iterations = 100,
          Bf_calculation = .x$Bf_calculation,
          tail = .x$tail,
          true_effect = .x$true_effect
        )
      }
    }
  ))

# Count occurrences of 1, 2, and 0 in hit column if output is not null
simulation_outputs <-
  simulation_outputs |>
  mutate(
    count_h1 = map_int(output, ~ {if (is.null(.x$error)) sum(.x$result$hit == 2) else NA_integer_}),
    count_h0 = map_int(output, ~ {if (is.null(.x$error)) sum(.x$result$hit == 1) else NA_integer_}),
    count_inconclusive = map_int(output, ~ {if (is.null(.x$error)) sum(.x$result$hit == 0) else NA_integer_})
  )

usethis::use_data(simulation_outputs, overwrite = TRUE)

simulation_results <-
  simulation_outputs |>
  unnest(cols = input) |>
  select(sd_of_theory, sd1, tpr, threshold, stopping_rule, Bf_calculation, tail, true_effect, n, hypothesis, count_h1, count_h0, count_inconclusive)

write_csv(simulation_results, here::here("inst/extdata/simulation_results.csv"))
