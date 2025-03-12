# Library ----------------------------------------------------------------------
library(tidyverse)
library(future.apply)
library(here)
devtools::load_all()

# Functions --------------------------------------------------------------------
# Create a function to safely run the simulation
safe_simulate_Bf <- purrr::safely(simulate_Bf)
Bf_samplesize_safe <- purrr::safely(Bf_samplesize)

# Run calculations -------------------------------------------------------------
# Calculation configurations ---------------------------------------------------
# Dataframe of the input parameters for the simulation
simulation_options <-
  expand.grid(
    sd_of_theory = 30,
    sd1 = 87.34,
    tpr = c(.5, .8, .9, .95),
    threshold = c(3, 6, 10),
    # stopping_rule = c("fixed", "optional"),
    stopping_rule = "fixed",
    # bf_calculation = c("normal", "cauchy"),
    bf_calculation = "normal",
    # tail = c(1, 2),
    tail = 1
    # true_effect = c(TRUE, FALSE)
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
      if (.x$bf_calculation == "normal") {
        Bf_samplesize_safe(
          distribution_prior = "normal",
          distribution_likelihood = "normal",
          sd_of_theory = 30,
          sd1 = 87.34,
          sd2 = NULL,
          tail = .x$tail,
          threshold = .x$threshold,
          tpr = .x$tpr,
          n_start = 5,
          n_end = 500,
          step = 1
        )
      } else if (.x$bf_calculation == "cauchy") {
        Bf_samplesize_safe(
          distribution_prior = "cauchy",
          distribution_likelihood = "normal",
          sd_of_theory = 30,
          sd1 = 87.34,
          sd2 = NULL,
          tail = .x$tail,
          threshold = .x$threshold,
          tpr = .x$tpr,
          n_start = 5,
          n_end = 500,
          step = 1
        )
      } else {
        stop("Unknown bf_calculation value")
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

simulation_samplesize_df <- 
  simulation_samplesize |> 
  unnest_wider(input) |> 
  select(-samplesize)

writexl::write_xlsx(simulation_samplesize_df, "inst/extdata/simulation_samplesize_bf.xlsx")

usethis::use_data(simulation_samplesize, overwrite = TRUE)

# Run simulations
simulation_outputs <-
  simulation_samplesize |>
  mutate(output = pmap(
    list(input, n, hypothesis), 
    ~ {
      if (is.na(..2)) {
        NULL
      } else {
        print(paste("Running simulation:", "sd_of_theory:", ..1$sd_of_theory, "sd:", ..1$sd1, "threshold:", ..1$threshold,
                    "stopping_rule:", ..1$stopping_rule, "n:", ..2, "bf_calculation:", ..1$bf_calculation,
                    "tail:", ..1$tail, "true_effect:", ifelse(..3 == "h1", TRUE, FALSE), "hypothesis:", ..3))
        
        safe_simulate_Bf(
          sd_of_theory = ..1$sd_of_theory,
          sd1 = ..1$sd1,
          sd2 = NULL,
          correlation = 0,
          threshold = ..1$threshold,
          stopping_rule = ..1$stopping_rule,
          n = ..2, 
          iterations = 100,
          bf_calculation = ..1$bf_calculation,
          tail = ..1$tail,
          true_effect = ifelse(..3 == "h1", TRUE, FALSE)
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
    count_inconclusive = map_int(output, ~ {if (is.null(.x$error)) sum(.x$result$hit == 0) else NA_integer_}),
    bf = map(output, ~{if (is.null(.x$error)) .x$result$bf else NA_real_})
  )

usethis::use_data(simulation_outputs, overwrite = TRUE)

simulation_results <-
  simulation_outputs |>
  unnest(cols = input) |>
  mutate(
    true_effect = ifelse(hypothesis == "h1", TRUE, FALSE)
  ) |> 
  select(sd_of_theory, sd1, tpr, threshold, stopping_rule, bf_calculation, tail, n, hypothesis, true_effect, count_h1, count_h0, count_inconclusive, bf)

write_csv(simulation_results, here::here("inst/extdata/simulation_results.csv"))
