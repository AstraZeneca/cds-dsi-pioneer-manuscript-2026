library(magrittr)
library(tidyverse)
library(rlang)
library(targets)
library(crew)
library(here)

source(here("r", "util.R"))
source(here("r", "prepare_analysis_data.R"))

source(here("r", "entimice", "entimice_functions.R"))

source(here("r", "breast", "util.R"))
source(here("r", "breast", "crcr.R"))
source(here("r", "breast", "prepare_analysis_data.R"))
source(here("r", "breast", "priors.R"))

tmp_dir <- file.path(Sys.getenv("TMPDIR"), "pioneer")  

tar_option_set(
  packages = c("tidyverse", "rlang", "here", "targets"),
  controller = crew_controller_local(workers = 12, seconds_timeout = 60, launch_max = 20),
  memory = "transient", garbage_collection = TRUE, 
  format = "qs",
  error = "continue"
)

lst(
  tar_target(util_stan_file, here("stan", "extern_util.stan"), format = "file"),
  tar_target(pfs_functions_file, here("stan", "extern_pfs_functions.stan"), format = "file"),
  tar_target(pfs_functions, cmdstan_expose_pfs_functions(util_stan_file, pfs_functions_file), deployment = "main", memory = "persistent", cue = tar_cue("always")),
 
  # Breast Data 
  
  tar_target(tumor_priors, get_tumor_priors()), 
  tar_target(pfs_priors, get_pfs_priors()),
  
  tar_target(db_data_path, "/wscratch/kmjq089/workspaces/dpo0083"),
  
  tar_target(# Breast clinical trial meta-data
    db_data_details, 
    tribble(
      ~study,          ~ idap,           ~ trial, 
      "d9673c00002",   "idap_20231204",  "Breast01",
      "d967hc00001",   "idap_20231129",  "Breast02", 
      "d9674c00001",   "idap_20231204",  "Breast03",
      "d9675c00001",   "idap_20231124",  "Breast04"
    ) |> 
      mutate(across(everything(), factor))),
  
  tar_target(db_sdtm_dm, read_all_trial_entimice_data(db_data_details, "deid_dm", "sdtm", team_dir = db_data_path)), # Demographics
  tar_target(db_sdtm_tu, read_all_trial_entimice_data(db_data_details, "deid_tu", "sdtm", team_dir = db_data_path)),
  tar_target(db_adam_adtr, read_all_trial_adtr_entimice_data(db_data_details, db_sdtm_tu, team_dir = db_data_path)), # Tumors
  tar_target(db_adam_adsl, read_all_trial_entimice_data(db_data_details, "deid_adsl", "adam", team_dir = db_data_path)),
  tar_target(db_adam_adcm, read_all_trial_entimice_data(db_data_details, "deid_adcm", "adam", team_dir = db_data_path)),
  tar_target(db_adam_adtte, read_all_trial_adtte_entimice_data(db_data_details, team_dir = db_data_path)),
  tar_target(db_adam_adrs, read_all_trial_entimice_data(db_data_details, "deid_adrs", "adam", team_dir = db_data_path)),
  
  tar_target(db_adam_adeg, read_all_trial_entimice_data(db_data_details, "deid_adeg", "adam", team_dir = db_data_path)), # ECG
  tar_target(db_adam_adis, read_all_trial_entimice_data(db_data_details, "deid_adis", "adam", team_dir = db_data_path)), # Immunogenicity/antibodies
  tar_target(db_adam_adlb, read_all_trial_entimice_data(db_data_details, "deid_adlb", "adam", team_dir = db_data_path)), # Lab analysis
  tar_target(db_adam_admh, read_all_trial_entimice_data(db_data_details, "deid_admh", "adam", team_dir = db_data_path)), # Medical history
  tar_target(db_adam_advs, read_all_trial_entimice_data(db_data_details, "deid_advs", "adam", team_dir = db_data_path)), # Vital signs
  
  tar_target(db_tumor_analysis_data, prepare_analysis_data(db_adam_adtr, db_adam_adtte, db_adam_adsl, db_adam_adcm, db_adam_adrs, pfs_functions), deployment = "main"),
  tar_target(unfiltered_db_pfs_analysis_data, prepare_raw_pfs_analysis_data(db_tumor_analysis_data, pfs_functions), deployment = "main"),
  tar_target(db_pfs_analysis_data, filter_pfs_analysis_data(unfiltered_db_pfs_analysis_data)), 
  tar_target(db_orr_analysis_data, prepare_orr_analysis_data(db_tumor_analysis_data)), 
  tar_target(treated_db_pfs_analysis_data, filter(db_pfs_analysis_data, treated)),
  tar_target(treated_db_pfs_stan_data, prepare_pfs_stan_data(treated_db_pfs_analysis_data, tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  tar_target(treated_trial_db_pfs_stan_data, 
             treated_db_pfs_analysis_data |> 
               nest(.by = trial) |> 
               deframe() |> 
               imap(\(d, tr) prepare_pfs_stan_data(mutate(d, trial = tr), tumor_priors, pfs_priors, pfs_functions)), deployment = "main"),
  tar_target(treated_db_obs_km, get_treated_obs_km(treated_trial_db_pfs_stan_data, pfs_functions), deployment = "main"),
  tar_target(db_early_tumor_pairs, get_early_tumor_pairs(treated_db_pfs_stan_data, pfs_functions), deployment = "main"),  
  tar_target(db_interval_censoring_data, get_ic_data(db_pfs_analysis_data, tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  
  tar_target(db_km_res_calculated, get_km_res(db_pfs_analysis_data, "pfs", tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  tar_target(db_km_res_from_data, get_km_res(db_pfs_analysis_data, "progress_week", tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  
  tar_target(db_confirmed_resp_analysis_data, prepare_confirmed_resp_analysis_data(db_tumor_analysis_data)),
  tar_target(treated_db_confirmed_resp_analysis_data, filter(db_confirmed_resp_analysis_data, treated)),
  # tar_target(
  #   treated_db_confirmed_resp_stan_data, 
  #   prepare_trial_confirmed_resp_stan_data(treated_db_confirmed_resp_analysis_data, confirmed_resp_priors, tumor_priors, pfs_conf_resp_priors, pfs_functions) |> 
  #     prepare_confirmed_resp_obs_km(pfs_functions),
  #   deployment = "main"
  # ),
  tar_target(db_km_confirmed_response, 
             with(treated_db_confirmed_resp_analysis_data, 
                  pfs_functions$estimate_kaplan_meier(confirmed_response_week, confirmed_response_censored, max(confirmed_response_week))),
             deployment = "main"),
  
  # tar_target(db_km_trial_confirmed_response, prepare_confirmed_resp_km(treated_confirmed_resp_stan_data), deployment = "main"),
  
  tar_target(db_km_confirmed_response_calendar, 
             with(treated_db_confirmed_resp_analysis_data, 
                  pfs_functions$estimate_kaplan_meier(
                    confirmed_response_week + experiment_start_week - 1, 
                    confirmed_response_censored, 
                    max(confirmed_response_week + experiment_start_week - 1))),
             deployment = "main"),
  
  tar_target(all_db_confirmed_resp_covar_formula, prepare_confirmed_resp_covar_formula()),
  tar_target(all_db_confirmed_resp_analysis_data, filter_conf_resp_missing_covar(treated_db_confirmed_resp_analysis_data, all_db_confirmed_resp_covar_formula)),
  tar_target(
    all_db_confirmed_resp_stan_data,
    prepare_confirmed_resp_stan_data(
      all_db_confirmed_resp_covar_formula, all_db_confirmed_resp_analysis_data, pfs_functions = pfs_functions
    ), 
    deployment = "main"
  ),
  
  tar_target(
    early_tumor_sums,
    with(
      all_db_confirmed_resp_stan_data, 
      pfs_functions$prepare_early_tumor_sums_covar(
        tumor_size, n_patient_tumors, n_measures, t_measure, 
        n_screening_t = pfs_functions$calc_n_screening_t(n_patient_tumors, n_measures, t_measure), 
        max_measures = 2
      )
    ),
    deployment = "main"
  ),
)
