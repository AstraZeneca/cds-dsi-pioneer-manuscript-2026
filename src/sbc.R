"
Usage: sbc.R [-c <cores>] [-n <num-sim>] [--append]

-c  Number of available cores [default: 12]
-n  Number of simulations to run [default: 12]
" |> 
  docopt::docopt(
    args = if (interactive()) "-c 3 -n 3" else commandArgs(TRUE),
  ) |> 
  purrr::modify_at(c("cores", "num_sim"), as.integer) -> cl_args

library(tidyverse)
library(rlang)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(here)

source(here("src", "util.R"))

options(mc.cores = cl_args$cores %/% 4)
future::plan(future::multisession(workers = cl_args$cores %/% 4))

fake_tumor_data <- read_rds(here("temp", "data", "fake_tumor.rds"))

pfs_model <- cmdstan_model(here("src", "pfs.stan"))

pfs_test_data <- lst(
    fit_data = FALSE,
    n_patients = 1000,
    gen_pfs = TRUE,
    max_pfs = 20,
    n_measures = 2,
    n_patient_tumors = count(fake_tumor_data, patient_id) |> pull(n),
    tumor_size = select(fake_tumor_data, tumor_size_1:tumor_size_2),
    pfs = rep(max_pfs, n_patients),
    censored = rep(1, n_patients)
  )

# Sample from the prior; no data.
pfs_res <- pfs_model$sample(data = pfs_test_data, refresh = 0)

sbc_data <- pfs_res |> # Get data from prior
  spread_rvars(rep_pfs[patient_index], rep_censored[patient_index]) |> 
  unnest_rvars() |> 
  ungroup() |> 
  filter(.draw <= cl_args$num_sim) |> 
  select(.draw, pfs = rep_pfs, censored = rep_censored) |> 
  nest(sim_data = !.draw) |>
  left_join( # Get the parameters that generated that data
    pfs_res |> 
      spread_rvars(tumor_stim_intercept, tumor_stim_coef[t]) |> 
      pivot_wider(names_from = t, values_from = tumor_stim_coef, names_prefix = "tumor_stim_coef_") |> 
      unnest_rvars() |> 
      ungroup() |> 
      select(!c(.iteration, .chain)),
    by = ".draw"
  ) |>
  mutate(
    furrr::future_map_dfr(.progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
      sim_data,
      # We need thinning when doing SBC using MCMC to break the correlation between samples.
      \(d) fit_sim_data(d, thin = 4) |> # Get posterior draws from simulation fit 
        spread_rvars(tumor_stim_intercept, tumor_stim_coef, ndraws = 1000) |> 
        rename_with(\(n) str_c("est_", n))
    ),
   
    # Rank statistics 
    r_tumor_stim_intercept = sum(est_tumor_stim_intercept < tumor_stim_intercept),
    r_tumor_stim_coef_1 = c(sum(est_tumor_stim_coef[, 1] < tumor_stim_coef_1)),
    r_tumor_stim_coef_2 = c(sum(est_tumor_stim_coef[, 2] < tumor_stim_coef_2)),
  )

if (cl_args$append) {
  try(sbc_data <- bind_rows(read_rds(here("temp", "data", "sbc.rds")), sbc_data))
} 

write_rds(sbc_data, here("temp", "data", "sbc.rds"))
