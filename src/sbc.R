"
Usage: sbc.R [-c <cores>] [-n <num-sim>]

-c  Number of available cores [default: 3]
-n  Number of simulations to run [default: 3]
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

options(mc.cores = cl_args$cores)
future::plan(future::multisession(workers = cl_args$cores))

fake_tumor_data <- read_rds(here("temp", "data", "fake_tumor.rds"))

pfs_model <- cmdstan_model(here("src", "pfs.stan"))

pfs_test_data <- lst(
    fit_data = FALSE,
    early_tumors_only = TRUE,
    n_patients = 1000,
    gen_pfs = TRUE,
    n_measures = 2,
    n_patient_tumors = count(fake_tumor_data, patient_id) |> pull(n),
    tumor_size = select(fake_tumor_data, tumor_size_1:tumor_size_2),
    pfs = rep(max_pfs, n_patients),
    right_censored = rep(1, n_patients),
    interval_censored = rep(0, n_patients)
  )

pfs_res <- pfs_model$sample(data = pfs_test_data, refresh = 0)

sbc_data <- pfs_res |> 
  spread_rvars(rep_pfs[patient_index], rep_right_censored[patient_index]) |> 
  unnest_rvars() |> 
  ungroup() |> 
  filter(.draw <= cl_args$num_sim) |> 
  select(.draw, pfs = rep_pfs, right_censored = rep_right_censored) |> 
  mutate(interval_censored = 0) |> 
  nest(sim_data = !.draw) |>
  left_join(
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
      \(d) fit_sim_data(d, thin = 4) |> 
        spread_rvars(tumor_stim_intercept, tumor_stim_coef, ndraws = 1000) |> 
        rename_with(\(n) str_c("est_", n))
    ),
    
    r_tumor_stim_intercept = sum(est_tumor_stim_intercept < tumor_stim_intercept),
    r_tumor_stim_coef_1 = c(sum(est_tumor_stim_coef[, 1] < tumor_stim_coef_1)),
    r_tumor_stim_coef_2 = c(sum(est_tumor_stim_coef[, 2] < tumor_stim_coef_2)),
  )

write_rds(sbc_data, here("temp", "data", "sbc.rds"))