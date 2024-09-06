library(conflicted)
library(magrittr)
library(tidyverse)
library(rlang)
library(targets)
library(stantargets)
library(crew)
library(here)
library(cmdstanr)

source(here("r", "util.R"))
source(here("r", "posterior.R"))
source(here("r", "crcr.R"))
source(here("r", "entimice_functions.R"))
source(here("r", "prepare_analysis_data.R"))
source(here("r", "priors.R"))

tmp_dir <- file.path(Sys.getenv("TMPDIR"), "adc-early-predict") |> 
  str_replace("^/scratch", "/wscratch") # Some SLURM nodes use the old TMPDIR /scratch

tar_option_set(
  packages = c("conflicted", "tidyverse", "rlang", "cmdstanr", "tidybayes", "here", "targets", "posterior"),
  controller = crew_controller_local(workers = 12, seconds_timeout = 60, launch_max = 20),
  memory = "transient", garbage_collection = TRUE, 
  format = "qs",
  error = "continue"
)
tar_config_set(store = file.path(tmp_dir, "_targets"))

options(cmdstanr_output_dir = file.path(tmp_dir, "fit"))

lst(
  # Prior hyperparameters
  
  tar_target(tumor_priors, get_tumor_priors()), 
  tar_target(pfs_priors, get_pfs_priors()),
  tar_target(confirmed_resp_priors, get_confirmed_resp_priors()),
  tar_target(pfs_conf_resp_priors, get_pfs_conf_resp_priors()),
  
  # Stan models
 
  tar_target(tumor_model_file, here("stan", "tumor.stan"), format = "file"),
  tar_target(pfs_model_file, here("stan", "pfs.stan"), format = "file"),
  tar_target(pfs2_model_file, here("stan", "pfs2.stan"), format = "file"),
  tar_target(pfs_orr_model_file, here("stan", "pfs_orr.stan"), format = "file"),
  tar_target(pfs_cr_model_file, here("stan", "pfs-confirmed-response.stan"), format = "file"),
  tar_target(crcr_model_file, here("stan", "crcr", "confresp-comprisk.stan"), format = "file"),
  tar_target(util_stan_file, here("stan", "extern_util.stan"), format = "file"),
  tar_target(pfs_functions_file, here("stan", "extern_pfs_functions.stan"), format = "file"),
  tar_target(recruit_model_file, here("stan", "recruit", "recruit.stan"), format = "file"),
  tar_target(recruit_maturity_model_file, here("stan", "recruit", "recruit_sample_maturity.stan"), format = "file"),
  tar_target(pfs_model, cmdstan_model(pfs_model_file)), 
  tar_target(pfs2_model, cmdstan_model(pfs2_model_file)), 
  # This builds the functions in Stan and makes them available in R. Now if you need these functions in downstream targets, you need to set 
  # cue to be "always"; loading a saved pfs_functions object will not work. Also, it needs to be on the main process as the targets that use it, hence, 
  # we use deployment = "main".
  tar_target(pfs_functions, cmdstan_expose_pfs_functions(util_stan_file, pfs_functions_file), deployment = "main", memory = "persistent"), #cue = tar_cue("always")),
  tar_target(tumor_model, cmdstan_model(tumor_model_file)),
  tar_target(pfs_orr_model, cmdstan_model(pfs_orr_model_file)),
  tar_target(pfs_cr_model, cmdstan_model(pfs_cr_model_file, force_recompile = TRUE, cpp_options = list(stan_threads = TRUE), stanc_options = list("O1"))),
  tar_target(crcr_model, cmdstan_model(crcr_model_file, force_recompile = TRUE, cpp_options = list(stan_threads = TRUE), stanc_options = list("O1"))),
  tar_target(recruit_model, cmdstan_model(recruit_model_file, force_recompile = TRUE)),
  tar_target(recruit_maturity_model, cmdstan_model(recruit_maturity_model_file, force_recompile = TRUE)),
 
  # Data 
  
  tar_target(# Clinical trial meta-data
    data_details, 
    tribble(
      ~study,          ~ idap,           ~ trial, 
      "d9673c00002",   "idap_20231204",  "Breast01",
      "d967hc00001",   "idap_20231129",  "Breast02", 
      "d9674c00001",   "idap_20231204",  "Breast03",
      "d9675c00001",   "idap_20231124",  "Breast04"
    ) |> 
      mutate(across(everything(), factor))),
  
  tar_target(sdtm_dm, read_all_trial_entimice_data(data_details, "deid_dm", "sdtm")), # Demographics
  tar_target(sdtm_tu, read_all_trial_entimice_data(data_details, "deid_tu", "sdtm")),
  tar_target(adam_adtr, read_all_trial_adtr_entimice_data(data_details, sdtm_tu)), # Tumors
  tar_target(adam_adsl, read_all_trial_entimice_data(data_details, "deid_adsl", "adam")),
  tar_target(adam_adcm, read_all_trial_entimice_data(data_details, "deid_adcm", "adam")),
  tar_target(adam_adtte, read_all_trial_adtte_entimice_data(data_details)),
  tar_target(adam_adrs, read_all_trial_entimice_data(data_details, "deid_adrs", "adam")),
  
  tar_target(adam_adeg, read_all_trial_entimice_data(data_details, "deid_adeg", "adam")), # ECG
  tar_target(adam_adis, read_all_trial_entimice_data(data_details, "deid_adis", "adam")), # Immunogenicity/antibodies
  tar_target(adam_adlb, read_all_trial_entimice_data(data_details, "deid_adlb", "adam")), # Lab analysis
  tar_target(adam_admh, read_all_trial_entimice_data(data_details, "deid_admh", "adam")), # Medical history
  tar_target(adam_advs, read_all_trial_entimice_data(data_details, "deid_advs", "adam")), # Vital signs
  
  tar_target(tumor_analysis_data, prepare_analysis_data(adam_adtr, adam_adtte, adam_adsl, adam_adcm, adam_adrs, pfs_functions), deployment = "main"),
  tar_target(unfiltered_pfs_analysis_data, prepare_raw_pfs_analysis_data(tumor_analysis_data, pfs_functions), deployment = "main"),
  tar_target(pfs_analysis_data, filter_pfs_analysis_data(unfiltered_pfs_analysis_data)), 
  tar_target(orr_analysis_data, prepare_orr_analysis_data(tumor_analysis_data)), 
  tar_target(treated_pfs_analysis_data, filter(pfs_analysis_data, treated)),
  tar_target(treated_pfs_stan_data, prepare_pfs_stan_data(treated_pfs_analysis_data, tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  tar_target(treated_trial_pfs_stan_data, 
             treated_pfs_analysis_data |> 
               nest(.by = trial) |> 
               deframe() |> 
               imap(\(d, tr) prepare_pfs_stan_data(mutate(d, trial = tr), tumor_priors, pfs_priors, pfs_functions)), deployment = "main"),
  tar_target(treated_obs_km, get_treated_obs_km(treated_trial_pfs_stan_data, pfs_functions), deployment = "main"),
  tar_target(early_tumor_pairs, get_early_tumor_pairs(treated_pfs_stan_data, pfs_functions), deployment = "main"),  
  tar_target(interval_censoring_data, get_ic_data(pfs_analysis_data, tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  
  tar_target(km_res_calculated, get_km_res(pfs_analysis_data, "pfs", tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  tar_target(km_res_from_data, get_km_res(pfs_analysis_data, "progress_week", tumor_priors, pfs_priors, pfs_functions), deployment = "main"),
  
  tar_target(recruit_phi, 7.5),
  
  tar_target(recruit_stan_data, treated_pfs_stan_data |> list_assign(n_trial_sim = 100)),
  tar_target(recruit_maturity_stan_data, 
             treated_confirmed_resp_stan_data |>
               rowwise() |> 
               mutate(stan_data = list(list_assign(stan_data, phi = recruit_phi, lambda = seq(25, 75, 10), pred_week = seq(6, 48, 6)) %>%
                                         list_assign(n_lambda = length(.$lambda), n_pred_week = length(.$pred_week))))),
  
  tar_target(confirmed_resp_analysis_data, prepare_confirmed_resp_analysis_data(tumor_analysis_data)),
  tar_target(treated_confirmed_resp_analysis_data, filter(confirmed_resp_analysis_data, treated)),
  tar_target(
    treated_confirmed_resp_stan_data, 
    prepare_trial_confirmed_resp_stan_data(treated_confirmed_resp_analysis_data, confirmed_resp_priors, tumor_priors, pfs_conf_resp_priors, pfs_functions) |> 
      prepare_confirmed_resp_obs_km(pfs_functions),
    deployment = "main"
  ),
  tar_target(km_confirmed_response, 
             with(treated_confirmed_resp_analysis_data, 
                  pfs_functions$estimate_kaplan_meier(confirmed_response_week, confirmed_response_censored, max(confirmed_response_week))),
             deployment = "main"),
  
  tar_target(km_trial_confirmed_response, prepare_confirmed_resp_km(treated_confirmed_resp_stan_data), deployment = "main"),
  
  tar_target(km_confirmed_response_calendar, 
             with(treated_confirmed_resp_analysis_data, 
                  pfs_functions$estimate_kaplan_meier(
                    confirmed_response_week + experiment_start_week - 1, 
                    confirmed_response_censored, 
                    max(confirmed_response_week + experiment_start_week - 1))),
             deployment = "main"),
  
  tar_target(all_confirmed_resp_covar_formula, prepare_confirmed_resp_covar_formula()),
  tar_target(all_confirmed_resp_analysis_data, filter_conf_resp_missing_covar(treated_confirmed_resp_analysis_data, all_confirmed_resp_covar_formula)),
  tar_target(
    all_confirmed_resp_stan_data,
    prepare_confirmed_resp_stan_data(
      all_confirmed_resp_covar_formula, all_confirmed_resp_analysis_data, 
      confirmed_resp_priors, .tumor_priors = tumor_priors, .pfs_priors = pfs_conf_resp_priors, 
      pfs_functions = pfs_functions
    ), 
    deployment = "main"
  ),
  
  # Experiment recruitment 
  
  tar_target(
    prior_recruit_fit, 
    recruit_model$sample(recruit_stan_data |> list_assign(fit_data = FALSE), 
                         parallel_chains = 4, iter_warmup = 200, iter_sampling = 200, 
                         output_basename = "prior_recruit", output_dir = file.path(tmp_dir, "fit")) |> 
      recover_types(select(treated_pfs_analysis_data, trial))
  ),
  
  tar_target(
    recruit_fit, 
    recruit_model$sample(recruit_stan_data,
                         parallel_chains = 4, iter_warmup = 200, iter_sampling = 200, 
                         output_basename = "recruit", output_dir = file.path(tmp_dir, "fit")) |> 
      recover_types(select(treated_pfs_analysis_data, trial))
  ),
  
  tar_target(recruit_maturity_res,
             rowwise(recruit_maturity_stan_data) |>
             mutate(fit = list(recruit_maturity_model$sample(
               stan_data,
               iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, # init = init_fun,
               output_dir = file.path(tmp_dir, "fit"), output_basename = str_c("recruit_maturity", trial)
             )
           )),
           pattern = map(recruit_maturity_stan_data)
  ),
  
  tar_target(recruit_maturity, get_sample_maturity_rvar(recruit_maturity_res)),
  
  # Tumor simulation 
   
  tar_target(n_patients, 250),
  tar_target(patient_measures, 41),
  tar_target(max_pfs, 45),
  tar_target(t_offset, 5),
  
  tar_target(tumor_sim_stan_data, get_sim_tumor_stan_data(n_patients, patient_measures, t_offset, tumor_priors)),
  tar_target(tumor_prior_res, tumor_model$sample(tumor_sim_stan_data, parallel_chains = 4, output_basename = "tumor_prior", output_dir = file.path(tmp_dir, "fit"))),
  tar_target(fake_tumor_data, get_fake_tumor_data(tumor_prior_res, patient_measures, t_offset, tumor_sim_stan_data)),
  tar_target(tumor_ppc_draws, get_tumor_ppc_draws(tumor_prior_res, patient_measures, t_offset, tumor_sim_stan_data), format = "fst_tbl"),
  tar_target(tumor_prior_fake_data, get_tumor_prior_fake_data(tumor_prior_res, patient_measures, t_offset, tumor_sim_stan_data)),
  tar_target(tumor_prior_fake_data_fit, 
             sim_tumor_fake_data(tumor_prior_fake_data, tumor_sim_stan_data, tumor_model, file.path(tmp_dir, "fit")), 
             pattern = map(tumor_prior_fake_data)),
  tar_target(tumor_prior_fake_data_rvar, get_tumor_prior_fake_data_rvar(tumor_prior_fake_data_fit, tumor_prior_res)),
  
  # Fitting tumor data
  
  tar_target(treated_tumor_analysis_data, filter(pfs_analysis_data, treated)),
  tar_target(tumor_stan_data, get_tumor_analysis_stan_data(treated_tumor_analysis_data, tumor_priors)),
  tar_target(n_full_measures, with(tumor_stan_data, n_measures + pfs_functions$calculate_n_missing_measures(n_measures, t_measure, n_patient_tumors)), deployment = "main"),
  tar_target(tumor_fit_res, 
             tumor_model$sample(tumor_stan_data, iter_warmup = 200, iter_sampling = 200, max_treedepth = 15, parallel_chains = 4, 
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "tumor")),
  tar_target(all_tumor_size, get_all_tumor_size(tumor_fit_res, treated_tumor_analysis_data, n_full_measures)),
  
  # PFS model prior prediction
  
  tar_target(early_fake_tumors, get_early_fake_tumors(fake_tumor_data)),
  tar_target(pfs_prior_stan_data, get_pfs_prior_stan_data(fake_tumor_data, early_fake_tumors, tumor_sim_stan_data, n_patients, max_pfs, pfs_priors)), 
  tar_target(max_t, get_max_t(pfs_prior_stan_data)),
  tar_target(pfs_prior_res, 
             pfs_model$sample(pfs_prior_stan_data, iter_warmup = 500, iter_sampling = 500, parallel_chains = 4, 
                              output_dir = file.path(tmp_dir, "fit"), output_basename = "pfs_prior")), 
  
  tar_target(rep_pfs_rvar, spread_rvars(pfs_prior_res, rep_pfs[patient_index], rep_right_censored[patient_index])),
  tar_target(rep_pfs_param_rvar, get_pfs_params(pfs_prior_res)), 
  tar_target(lambda_rvar, get_lambda_rvar(pfs_prior_res)),
  tar_target(expected_pfs_rvar, get_expected_pfs_rvar(pfs_prior_res)),
  tar_target(prior_dp_prob_draws, get_disease_progress_prob(pfs_prior_res, pfs_prior_stan_data, max_t, n = 2), format = "fst_tbl"),
  
  # PFS fake data simulation
  
  tar_target(t_to_keep, c(-t_offset, 0, seq(1 - t_offset, patient_measures - t_offset - 1, by = 5)) |> sort()),
  tar_target(pfs_fake_data_sim_settings, list_assign(pfs_prior_stan_data, add_trial_level = FALSE)),
  tar_target(pfs_fake_data, get_pfs_sim_seed_draws(pfs_fake_data_sim_settings, 12, pfs_model)),
  tar_target(no_ic_lstm, with(pfs_fake_data_sim_settings, list_measures(t_measure, n_measures, n_patient_tumors))),
  tar_target(
    no_ic_fake_data_sim,
    run_simulation(pfs_fake_data$.draw, pfs_fake_data_sim_settings, pfs_fake_data$sim_data[[1]], no_ic_lstm, 
                   file.path(tmp_dir, "fit"), "fake_data_sim", 
                   FALSE, TRUE, pfs_model, keep_fit = TRUE), 
    pattern = map(pfs_fake_data)
  ), 
  tar_target(censored_pfs_fake_data_sim_settings, drop_missing_measures(pfs_fake_data_sim_settings, t_to_keep, TRUE)),
  tar_target(ic_lstm, with(censored_pfs_fake_data_sim_settings, list_measures(t_measure, n_measures, n_patient_tumors))),
  tar_target(
    ic_fake_data_sim,
    run_simulation(pfs_fake_data$.draw, censored_pfs_fake_data_sim_settings, pfs_fake_data$sim_data[[1]], ic_lstm, 
                   file.path(tmp_dir, "fit"), "fake_data_sim_ic", 
                   FALSE, TRUE, pfs_model, keep_fit = TRUE), 
    pattern = map(pfs_fake_data)
  ), 
  tar_target(
    ic_ignored_fake_data_sim,
    run_simulation(pfs_fake_data$.draw, censored_pfs_fake_data_sim_settings, pfs_fake_data$sim_data[[1]], ic_lstm, 
                   file.path(tmp_dir, "fit"), "fake_data_sim_ic_ignored", 
                   TRUE, TRUE, pfs_model, keep_fit = TRUE), 
    pattern = map(pfs_fake_data)
  ), 
  
  # TODO migrate SBC to targets
  
  # Fitting PFS with no covariates 
  
  tar_target(no_tumor_treated_pfs_stan_data, 
             list_assign(treated_pfs_stan_data, tumor_hazard_type = 0, add_trial_level = TRUE)),
  tar_target(prior_no_tumor_pfs_res, 
             no_tumor_treated_pfs_stan_data |> 
               list_modify(fit_data = FALSE) |>
               pfs_model$sample(parallel_chains = 4, iter_warmup = 300, iter_sampling = 300, refresh = 0, 
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_no_tumor_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  tar_target(no_tumor_pfs_res,
             no_tumor_treated_pfs_stan_data |> 
               pfs_model$sample(parallel_chains = 4, iter_warmup = 300, iter_sampling = 300, refresh = 0,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "no_tumor_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  tar_target(ic_ignored_no_tumor_pfs_res,
             no_tumor_treated_pfs_stan_data |>
              list_modify(ignore_interval_censoring = TRUE) |>
              pfs_model$sample(parallel_chains = 4, iter_warmup = 300, iter_sampling = 300, refresh = 0,
                               output_dir = file.path(tmp_dir, "fit"), output_basename = "ic_ignored_no_tumor_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  
  # Fitting PFS predicted by tumors count only 
  
  tar_target(tumor_count_treated_pfs_stan_data, 
             list_assign(treated_pfs_stan_data, tumor_hazard_type = 4, add_trial_level = TRUE)),
  tar_target(tumor_count_init, create_pfs_initializer(tumor_count_treated_pfs_stan_data)),
  
  tar_target(prior_tumor_count_pfs_res,
             tumor_count_treated_pfs_stan_data |> 
               list_assign(fit_data = FALSE, add_trial_level = TRUE) |>
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = tumor_count_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_tumor_count_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  tar_target(tumor_count_pfs_res,
             tumor_count_treated_pfs_stan_data |> 
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = tumor_count_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "tumor_count_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  
  # Fitting PFS predicted by first two tumor measures
  
  tar_target(two_tumor_treated_pfs_stan_data, 
             list_assign(treated_pfs_stan_data, tumor_hazard_type = 1, add_trial_level = TRUE, add_tumor_location_level = FALSE)),
  tar_target(two_tumor_init, create_pfs_initializer(two_tumor_treated_pfs_stan_data)),
  tar_target(prior_two_tumor_pfs_res,
             two_tumor_treated_pfs_stan_data |> 
               list_assign(fit_data = FALSE) |>
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = two_tumor_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_two_tumor_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  tar_target(two_tumor_pfs_res,
             two_tumor_treated_pfs_stan_data |> 
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = two_tumor_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "two_tumor_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  
  # Looking at partial pooling for PFS prediction
  
  tar_target(two_tumor_treated_trial_pfs_stan_data, 
             map(treated_trial_pfs_stan_data, 
                 \(s) list_assign(s, tumor_hazard_type = 2) |> modify_at(c("patient", "patient_trial"), factor)),
             iteration = "list"),
  tar_target(trial_names, names(two_tumor_treated_trial_pfs_stan_data)),
  tar_target(two_tumor_trial_pfs_res_list, 
             fit_by_trial(trial_names, two_tumor_treated_trial_pfs_stan_data, pfs_model), 
             pattern = map(two_tumor_treated_trial_pfs_stan_data, trial_names), iteration = "list"),
  tar_target(two_tumor_trial_pfs_res, set_names(two_tumor_trial_pfs_res_list, trial_names)),
  
  # Fitting PFS predicted by tumor size change
  
  tar_target(tumor_change_treated_pfs_stan_data, 
             list_modify(treated_pfs_stan_data, tumor_hazard_type = 3, add_trial_level = TRUE)),
  tar_target(tumor_change_init, create_pfs_initializer(tumor_change_treated_pfs_stan_data)),
  tar_target(prior_tumor_change_pfs_res, 
             tumor_change_treated_pfs_stan_data |> 
               list_modify(fit_data = FALSE) |>
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = tumor_change_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_tumor_change_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  tar_target(tumor_change_pfs_res,
             tumor_change_treated_pfs_stan_data |> 
               pfs_model$sample(iter_warmup = 300, iter_sampling = 300, init = tumor_change_init, parallel_chains = 4,
                                output_dir = file.path(tmp_dir, "fit"), output_basename = "tumor_change_pfs") |> 
               recover_types(treated_pfs_analysis_data)),
  
  # tar_target(tumor_design_matrix, get_tumor_design_matrix(treated_pfs_stan_data, pfs_functions)),
 
  # Model for confirmed response, pooled
  
  tar_target(cif_hb, seq(0, 1, by = 0.05)),
  tar_target(cif_t, seq(1, 100, 10)),
  tar_target(conf_resp_hb, seq(0, 144, by = 6)),
  tar_target(crcr_hazard_ratio_hb, seq(0, 3, by = 0.1)),
  
  tar_target(
    prior_confirmed_resp_comp_risk_res,
    crcr_model$sample(
      all_confirmed_resp_stan_data |> list_assign(fit_data = FALSE),
      iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, threads_per_chain = 4, refresh = 10, # init = init_fun,
      output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_confirmed_resp_comp_risk"
    )
  ),
  
  tar_target(
    confirmed_resp_comp_risk_res,
    crcr_model$sample(
      all_confirmed_resp_stan_data,
      iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, threads_per_chain = 4, refresh = 10, # init = init_fun,
      output_dir = file.path(tmp_dir, "fit"), output_basename = "confirmed_resp_comp_risk"
    )
  ),
   
  # Predicting survival using confirmed response and tumor sizes, pooled 
  
  tar_target(
    prior_all_confirmed_resp_pfs_res,
    pfs_cr_model$sample(
      all_confirmed_resp_stan_data |> list_assign(fit_data = FALSE),
      iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, threads_per_chain = 4, refresh = 5,
      output_dir = file.path(tmp_dir, "fit"), output_basename = "prior_confirmed_resp_pfs"
    )
  ),
  
  tar_target(
    all_confirmed_resp_pfs_res,
    pfs_cr_model$sample(
      all_confirmed_resp_stan_data,
      iter_warmup = 300, iter_sampling = 300, parallel_chains = 4, threads_per_chain = 4, refresh = 5,
      output_dir = file.path(tmp_dir, "fit"), output_basename = "confirmed_resp_pfs"
    )
  ),
  
  tar_target(prior_crcr_res, prior_all_confirmed_resp_pfs_res),
  tar_target(crcr_res, all_confirmed_resp_pfs_res),
  
  tar_target(
    prior_all_rep_confirmed_response_bindist, get_all_rep_confirmed_response_bindist(prior_crcr_res, all_confirmed_resp_stan_data, conf_resp_hb)
  ),
  tar_target(
    all_rep_confirmed_response_bindist, get_all_rep_confirmed_response_bindist(crcr_res, all_confirmed_resp_stan_data, conf_resp_hb)
  ),
  tar_target(
    prior_all_crcr_hazard_ratio_bindist, 
    get_all_conf_resp_hazard_ratios_bindist(prior_crcr_res, all_confirmed_resp_stan_data, crcr_hazard_ratio_hb)
  ),
  tar_target(
    all_crcr_hazard_ratio_bindist, 
    get_all_conf_resp_hazard_ratios_bindist(crcr_res, all_confirmed_resp_stan_data, crcr_hazard_ratio_hb)
  ),
  tar_target(
    prior_all_conf_resp_cif_bindist, 
    get_all_conf_resp_cif_bindist(prior_crcr_res, all_confirmed_resp_stan_data, cif_t, cif_hb, ndraws = 300)
  ),
  tar_target(
    all_conf_resp_cif_bindist, 
    get_all_conf_resp_cif_bindist(crcr_res, all_confirmed_resp_stan_data, cif_t, cif_hb, ndraws = 300)
  ),
  
  tar_target(
    prior_trial_conf_resp_median_pfs, 
    get_all_median_pfs_conf_resp(prior_all_confirmed_resp_pfs_res, all_confirmed_resp_analysis_data),
    storage = "worker", retrieval = "worker" 
  ),
  tar_target(
    trial_conf_resp_median_pfs, 
    get_all_median_pfs_conf_resp(all_confirmed_resp_pfs_res, all_confirmed_resp_analysis_data),
    storage = "worker", retrieval = "worker" 
  ),
  tar_target(prior_trial_conf_resp_km_est, get_all_pfs_conf_resp_km_est(prior_all_confirmed_resp_pfs_res, all_confirmed_resp_analysis_data)),
  tar_target(trial_conf_resp_km_est, get_all_pfs_conf_resp_km_est(all_confirmed_resp_pfs_res, all_confirmed_resp_analysis_data)),
  tar_target(prior_trial_pfs_conf_resp_hazard_ratio, get_all_pfs_conf_resp_log_hazard_ratio(prior_all_confirmed_resp_pfs_res)),
  tar_target(trial_pfs_conf_resp_hazard_ratio, get_all_pfs_conf_resp_log_hazard_ratio(all_confirmed_resp_pfs_res)),
  tar_target(prior_pfs_conf_resp_param, get_all_pfs_pred_param(prior_all_confirmed_resp_pfs_res), storage = "worker", retrieval = "worker"),
  tar_target(pfs_conf_resp_param, get_all_pfs_pred_param(all_confirmed_resp_pfs_res), storage = "worker", retrieval = "worker"),
  
  # Leave-one-trial-out bootstrap
  
  tar_target(leave_out_trials, unique(all_confirmed_resp_stan_data$patient_trial)),
  tar_target(
    bootstrap_settings,
    lst(
      n_bootstrap_sample = 30, 
      bootstrap_cr_maturity_rates = seq(0.05, 0.8, 0.05),
      n_bootstrap_cr_maturity_rates = length(bootstrap_cr_maturity_rates), 
      bootstrap_pfs_maturity_rates = bootstrap_cr_maturity_rates,
      n_bootstrap_pfs_maturity_rates = length(bootstrap_pfs_maturity_rates),
      n_fixed_bootstrap_samples = 10 
    )
  ),
  
  tar_target(
    bootstrap_confirmed_resp_pfs_res, 
    run_bootstrap_pfs_cr(pfs_cr_model, leave_out_trials, all_confirmed_resp_stan_data, bootstrap_settings), # iter_warmup = 2, iter_sampling = 2), 
    pattern = map(leave_out_trials)
  ),
  
  tar_target(
    pfs_conf_resp_bootstrap_cr_median_pfs, 
    get_pfs_conf_resp_bootstrap_cr_median_pfs(bootstrap_confirmed_resp_pfs_res), pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  tar_target(
    pfs_conf_resp_bootstrap_pfs_median_pfs, 
    get_pfs_conf_resp_bootstrap_pfs_median_pfs(bootstrap_confirmed_resp_pfs_res), pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  tar_target(
    fixed_bootstrap_cr_median_pfs, 
    get_fixed_bootstrap_cr_median_pfs(bootstrap_confirmed_resp_pfs_res), pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  tar_target(
    fixed_bootstrap_pfs_median_pfs, 
    get_fixed_bootstrap_pfs_median_pfs(bootstrap_confirmed_resp_pfs_res), pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  tar_target(
    pfs_conf_resp_bootstrap_param, 
    get_pfs_pred_param(bootstrap_confirmed_resp_pfs_res), pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  tar_target(
    bootstrap_cr_median_pfs_draws, 
    get_cr_median_pfs_draws(bootstrap_confirmed_resp_pfs_res, ndraws = 100),
    pattern = map(bootstrap_confirmed_resp_pfs_res),
    storage = "worker", retrieval = "worker" 
  ),
  
  # PFS CR separately
  
  tar_target(prior_confirmed_resp_pfs_res,
             treated_confirmed_resp_stan_data |> 
               rowwise() |> 
               mutate(
                 stan_data = stan_data |> list_assign(fit_data = FALSE) |> list(),
                 fit = list(pfs_cr_model$sample(
                   stan_data,
                   iter_warmup = 400, iter_sampling = 400, parallel_chains = 4, threads_per_chain = 4, init = init_fun,
                   output_dir = file.path(tmp_dir, "fit"), output_basename = str_c("prior_confirmed_resp_pfs_", trial)
                 )
             )), 
           pattern = map(treated_confirmed_resp_stan_data)
  ),
  
  tar_target(confirmed_resp_pfs_res,
             treated_confirmed_resp_stan_data |> 
               rowwise() |> 
               mutate(
                 fit = list(pfs_cr_model$sample(
                   stan_data,
                   iter_warmup = 400, iter_sampling = 400, parallel_chains = 4, threads_per_chain = 4, init = init_fun,
                   output_dir = file.path(tmp_dir, "fit"), output_basename = str_c("confirmed_resp_pfs_", trial)
                 )
             )), 
           pattern = map(treated_confirmed_resp_stan_data)
  ),
  
  # Extracting result rvars from model fit
  
  tar_target(prior_conf_resp_median_pfs, get_median_pfs_conf_resp(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  tar_target(conf_resp_median_pfs, get_median_pfs_conf_resp(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  tar_target(prior_conf_resp_km_est, get_pfs_conf_resp_km_est(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  tar_target(conf_resp_km_est, get_pfs_conf_resp_km_est(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  # tar_target(prior_pfs_conf_resp_marginal_exit_prob, 
  #            get_pfs_conf_resp_marginal_exit_prob(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  # tar_target(pfs_conf_resp_marginal_exit_prob, 
  #            get_pfs_conf_resp_marginal_exit_prob(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  # tar_target(prior_crcr_hazard_ratio, get_conf_resp_hazard_ratios(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  # tar_target(prior_rep_confirmed_response, get_rep_confirmed_response(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  # tar_target(rep_confirmed_response, get_rep_confirmed_response(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  # tar_target(crcr_hazard_ratio, get_conf_resp_hazard_ratios(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  tar_target(prior_pfs_conf_resp_hazard_ratio, get_pfs_conf_resp_log_hazard_ratio(prior_confirmed_resp_pfs_res), pattern = map(prior_confirmed_resp_pfs_res)),
  tar_target(pfs_conf_resp_hazard_ratio, get_pfs_conf_resp_log_hazard_ratio(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)),
  # tar_target(pfs_conf_resp_bootstrap_median_pfs, get_pfs_conf_resp_bootstrap_median_pfs(confirmed_resp_pfs_res), pattern = map(confirmed_resp_pfs_res)), 
  # tar_target(prior_conf_resp_cif, get_conf_resp_cif(prior_confirmed_resp_pfs_res, cif_t, cif_hb), pattern = map(prior_confirmed_resp_pfs_res)),
  # tar_target(conf_resp_cif, get_conf_resp_cif(confirmed_resp_pfs_res, cif_t, cif_hb), pattern = map(confirmed_resp_pfs_res)),
)
