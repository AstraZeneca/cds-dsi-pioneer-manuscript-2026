"
Usage: sbc.R <cores> <num-sim> <output-name> [--append] [--censor-intervals=<intervals> --keep-only --ignore-interval-censoring]

--censor-intervals=<intervals>  Intervals to censor in the data
--keep-only  Keep only the intervals provided for censoring
" |> 
  docopt::docopt(
    args = if (interactive()) "12 3 test  --censor-intervals=0,6,12,18,24,30,36 --keep-only" else commandArgs(TRUE),
    # args = if (interactive()) "12 12 test" else commandArgs(TRUE),
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

max_pfs <- 45 
tmp_dir <- file.path(Sys.getenv("TMPDIR"), "adc-early-predict") |> 
  str_replace("^/scratch", "/wscratch") # Some SLURM nodes use the old TMPDIR /scratch

cat("Temporary folder:", tmp_dir, "\n")

source(here("r", "util.R"))
source(here("r", "priors.R"))

future::plan(future::multisession(workers = cl_args$cores %/% 4))

# Load pre-generated tumor data
tumor_test_data <- rjson::fromJSON(file = file.path(tmp_dir, "data", "prior_tumor.json"))
fake_tumor_data <- read_rds(file.path(tmp_dir, "data", "fake_tumor.rds"))

pfs_model <- cmdstan_model(here("stan", "pfs.stan"))

# A prior only run to generate datasets
pfs_test_data <- tumor_test_data |> 
  list_modify(
    fit_data = FALSE,
    use_tumor_model = FALSE,
    tumor_hazard_type = 1,
    ignore_interval_censoring = FALSE,
    gen_pfs = TRUE,
    gen_interval_censored = TRUE,
    add_trial_level = FALSE, 
    add_tumor_location_level = FALSE, 
    fit_post_2nd_meaure_only = FALSE,
    
    grid_tumors = array(NA, dim = 0),
    n_grid_tumors = 0,
    tumor_size = fake_tumor_data$tumor_size,
    pfs = rep(max_pfs, tumor_test_data$n_patients),
    right_censored = rep(FALSE, tumor_test_data$n_patients),
    death_week = rep(0, tumor_test_data$n_patients),
    
    !!!tumor_priors,
    !!!pfs_priors,
  ) |> 
  drop_missing_measures(cl_args$censor_intervals, keep_only = cl_args$keep_only)  

sbc_data <- run_sbc_sims(pfs_model, pfs_test_data, cl_args$num_sim, ignore_interval_censoring = cl_args$ignore_interval_censoring, output_name = cl_args$output_name) |> 
  transmute(
    .draw,
   
    # Rank statistics
    # For each simulation calculate the number of posterior parameter samples that are less than the true parameter value. 
    r_tumor_stim_intercept = sum(est$tumor_stim_pop_intercept < true$tumor_stim_pop_intercept),
    r_tumor_stim_coef_1 = sum(est$tumor_stim_pop_coef_1 < true$tumor_stim_pop_coef_1),
    r_tumor_stim_coef_2 = sum(est$tumor_stim_pop_coef_2 < true$tumor_stim_pop_coef_2),
    r_log_lambda_gp_intercept = sum(est$log_lambda_gp_intercept < true$log_lambda_gp_intercept),
    r_log_lambda_gp_alpha = sum(est$log_lambda_gp_alpha < true$log_lambda_gp_alpha),
    r_log_lambda_gp_rho = sum(est$log_lambda_gp_rho < true$log_lambda_gp_rho),
    r_base_cond_expected_pfs = sum(est$base_cond_expected_pfs < true$base_cond_expected_pfs),
    r_one_tumor_cond_expected_pfs = sum(est$one_tumor_cond_expected_pfs < true$one_tumor_cond_expected_pfs),
    r_base_cond_median_pfs = sum(est$base_cond_median_pfs < true$base_cond_median_pfs),
    r_one_tumor_cond_median_pfs = sum(est$one_tumor_cond_median_pfs < true$one_tumor_cond_median_pfs),
  )

if (cl_args$append) {
  sbc_data <- try(bind_rows(read_rds(file.path(tmp_dir, "data", str_glue("{cl_args$output_name}.rds"))), sbc_data))
} 

write_rds(sbc_data, file.path(tmp_dir, "data", str_glue("{cl_args$output_name}.rds")))
