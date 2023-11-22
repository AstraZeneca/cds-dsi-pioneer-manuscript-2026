"
Usage: sbc.R <cores> <num-sim> <output-name> [--append] [--censor-intervals=<intervals> --ignore-interval-censoring]

-c                              Number of available cores 
-n                              Number of simulations to run 
-o                              Name to use for SBC run files, etc. 
--censor-intervals=<intervals>  Intervals to censor in the data
" |> 
  docopt::docopt(
    args = if (interactive()) "12 3 test --censor-intervals=3,4,6,9" else commandArgs(TRUE),
  ) -> cl_args

library(tidyverse)
library(rlang)
library(cmdstanr)
library(posterior)
library(tidybayes)
library(here)

cl_args <- cl_args |> 
  purrr::modify_at("censor_intervals", \(si) str_split(si, fixed(",")) |> list_c() |> unique()) |> 
  purrr::modify_at(c("cores", "num_sim", "censor_intervals"), as.integer) -> cl_args

max_pfs <- 20
tmp_dir <- file.path(Sys.getenv("TMP"), "adc-early-predict")

source(here("src", "util.R"))

options(mc.cores = cl_args$cores %/% 4)
future::plan(future::multisession(workers = cl_args$cores %/% 4))

tumor_test_data <- rjson::fromJSON(file = file.path(tmp_dir, "data", "prior_tumor.json"))
fake_tumor_data <- read_rds(file.path(tmp_dir, "data", "fake_tumor.rds"))

pfs_model <- cmdstan_model(here("src", "pfs.stan"))

pfs_test_data <- tumor_test_data |> 
  list_modify(
    fit_data = FALSE,
    early_tumors_only = TRUE,
    ignore_interval_censoring = FALSE,
    gen_pfs = TRUE,
    gen_interval_censored = TRUE,
    tumor_size = fake_tumor_data$tumor_size,
    pfs = rep(max_pfs, tumor_test_data$n_patients),
    
    log_lambda_gp_intercept_mean = -3,
    log_lambda_gp_intercept_sd = 0.5,
    tumor_stim_intercept_sd = 0.5,
    tumor_stim_coef_sd = c(0.25, 0.25),
  ) |> 
  drop_missing_measures(cl_args$censor_intervals)  

# Sample from the prior: no data.
pfs_res <- pfs_model$sample(data = pfs_test_data, refresh = 0)

sbc_data <- pfs_res |> 
  spread_rvars(rep_pfs[patient_index]) |> 
  unnest_rvars() |> 
  ungroup() |> 
  filter(.draw <= cl_args$num_sim) |> 
  select(.draw, pfs = rep_pfs) |> 
  nest(sim_data = !.draw) |>
  left_join( # Get the parameters that generated that data
    pfs_res |> 
      spread_rvars(tumor_stim_intercept, tumor_stim_coef[t], log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho) |> 
      pivot_wider(names_from = t, values_from = tumor_stim_coef, names_prefix = "tumor_stim_coef_") |> 
      unnest_rvars() |> 
      ungroup() |> 
      select(!c(.iteration, .chain)),
    by = ".draw"
  ) |>
  mutate(
    furrr::future_map2_dfr(.draw, sim_data, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
      function(sim_id, d, cl_args) { 
        fit_sim_data(
          pfs_test_data,
          d, 
          gen_pfs = FALSE, 
          thin = 4, # We need thinning when doing SBC using MCMC to break the correlation between samples.
          output_basename = str_glue("{cl_args$output_name}_{sim_id}"),
          output_dir = file.path(tmp_dir, "fit"), 
          ignore_interval_censoring = cl_args$ignore_interval_censoring 
        ) |> # Get posterior draws from simulation fit 
          spread_rvars(tumor_stim_intercept, tumor_stim_coef,  log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho, ndraws = 1000) |> 
          rename_with(\(n) str_c("est_", n))
      }, 
      cl_args = cl_args
    ),
  ) |> 
  transmute(
    .draw,
   
    # Rank statistics 
    r_tumor_stim_intercept = sum(est_tumor_stim_intercept < tumor_stim_intercept),
    r_tumor_stim_coef_1 = c(sum(est_tumor_stim_coef[, 1] < tumor_stim_coef_1)),
    r_tumor_stim_coef_2 = c(sum(est_tumor_stim_coef[, 2] < tumor_stim_coef_2)),
    r_log_lambda_gp_intercept = sum(est_log_lambda_gp_intercept < log_lambda_gp_intercept),
    r_log_lambda_gp_alpha = sum(est_log_lambda_gp_alpha < log_lambda_gp_alpha),
    r_log_lambda_gp_rho = sum(est_log_lambda_gp_rho < log_lambda_gp_rho),
  )

if (cl_args$append) {
  sbc_data <- try(bind_rows(read_rds(file.path(tmp_dir, "data", str_glue("{cl_args$output_name}.rds"))), sbc_data))
} 

write_rds(sbc_data, file.path(tmp_dir, "data", str_glue("{cl_args$output_name}.rds")))
