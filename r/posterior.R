# nolint start: object_usage_linter

# PFS and Confirmed Response Analysis Functions
#
# This file contains a collection of functions for analyzing Progression-Free Survival (PFS)
# and Confirmed Response data in clinical trials. The functions cover various aspects of
# survival analysis, including:
#
# - Extracting and processing model results
# - Calculating Kaplan-Meier estimates
# - Handling bootstrap samples
# - Computing hazard ratios
# - Generating prediction parameters
# - Power scaling of variables
#
# These functions are designed to work with Stan model outputs and provide tools for
# comprehensive analysis of clinical trial data, focusing on PFS and confirmed response
# endpoints.

get_pfs_conf_resp_marginal_exit_prob <- function(res) {
  res |>
    ungroup() |> 
    transmute(trial, prob_rvars = map(fit, \(f) spread_rvars(f, marginal_exit_prob[i, k, t]))) |> 
    unnest(prob_rvars)
}

get_trial_sim_crcr_pfs <- function(res, stan_data) {
  spread_rvars(res, sim_pfs[i], sim_censored[i], forecast_pfs[i], forecast_censored[i]) |>
    mutate(usubjid = stan_data$patient)
}

get_sim_pfs_conf_resp <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      sim_pfs = map2(
        fit, analysis_data, 
        \(f, d) spread_rvars(f, sim_pfs[i], sim_censored[i]) |>
          left_join(transmute(d, i = seq(n()), pfs, right_censored), by = "i", relationship = "one-to-one")
      )
    ) |> 
    unnest(sim_pfs)
}

get_pfs_conf_resp_km_est <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |> 
    map_dfr(\(r) spread_rvars(r, km_est[t]), .id = "trial") 
}

get_all_pfs_conf_resp_km_est <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
  spread_rvars(res, trial_km_est[trial, t], forecast_trial_km_est[trial, t])
}

get_median_pfs_conf_resp <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, sim_median_pfs))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_median_pfs_conf_resp <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
   spread_rvars(res, sim_trial_median_pfs[trial], forecast_trial_median_pfs[trial]) 
}

get_pfs_n <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
   spread_rvars(res, sim_trial_pfs6[trial], sim_trial_pfs9[trial], forecast_trial_pfs6[trial], forecast_trial_pfs9[trial]) 
}

get_pfs_conf_resp_log_hazard_ratio <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, time_invariant_log_hazard_ratio[k, i]) |> 
                                  mutate(time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_pfs_conf_resp_hazard_ratio <- function(res, stan_data) {
  spread_rvars(res, time_invariant_log_hazard_ratio[k, i]) |>
    mutate(
      time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response")) 
    ) |> 
    left_join(as_tibble(stan_data["patient_trial"]) |> mutate(i = seq(n())), by = "i", relationship = "many-to-one") |> 
    rename(trial = patient_trial)  
}

get_all_pfs_conf_resp_hazard_ratio_bindist <- function(res, stan_data, hb) {
  get_all_pfs_conf_resp_hazard_ratio(res, stan_data) |> 
    group_by(trial, k) |> 
    reframe(t = hb[-length(hb)], bindist = rvar_sample_hist(time_invariant_hazard_ratio, hb))  
}

get_pfs_conf_resp_bootstrap_variables <- function(res, bs_sample1, bs_sample2, rates_name, ..., summarize = TRUE) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(
        spread_rvars(fit, ...) |>  
          mutate(
            n_bs_sample = {{ bs_sample1 }} + {{ bs_sample2 }},
            across(ends_with("predicted"), \(n)  n / n_bs_sample, .names = "{.col}_prop")
          ) %>% { 
            if (summarize) {
              unnest_rvars(.) |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
                point_interval(
                  na.rm = TRUE, # if n_bs_sample is 0, we'll get some NaNs. These are few so we'll bite the bullet and drop them.
                  .width = c(0.5, 0.8)
                ) 
            } else .
          } |> 
          left_join(as_tibble(stan_data[rates_name]) |> mutate(r = seq(n())), by = "r", relationship = "many-to-one")
      )
    ) |> 
    ungroup() |> 
    unnest(rv)  
}

get_pfs_conf_resp_bootstrap_cr_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, "bootstrap_cr_maturity_rates",
    bs_cr_prediction_calendar_week[r], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r],
    bs_cr_median_pfs[r], bs_cr_orr[r],
    n_bs_cr_conf_resp_predicted[r], n_bs_cr_pfs_predicted[r],
    summarize = FALSE
  ) |> 
    mutate(log_bs_cr_median_pfs = log(bs_cr_median_pfs)) |> 
    unnest_rvars() |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
    point_interval(na.rm = TRUE, .width = c(0.5, 0.8))
}

get_pfs_conf_resp_bootstrap_pfs_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, "bootstrap_pfs_maturity_rates",
    bs_cr_conf_resp_censored_prop[r], bs_pfs_prediction_calendar_week[r], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r],
    bs_pfs_median_pfs[r], bs_pfs_orr[r],
    n_bs_pfs_conf_resp_predicted[r], n_bs_pfs_pfs_predicted[r]
  )
}    

get_fixed_bootstrap_cr_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, "bootstrap_cr_maturity_rates",
    fixed_bs_cr_median_pfs[r, f], fixed_bs_cr_orr[r, f], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r],
    n_fixed_bs_cr_conf_resp_predicted[r, f], n_fixed_bs_cr_pfs_predicted[r, f],
    summarize = FALSE
  ) |> 
    select(!c(n_bs_sample, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified))
}

get_fixed_bootstrap_pfs_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, "bootstrap_pfs_maturity_rates",
    fixed_bs_pfs_median_pfs[r, f], fixed_bs_pfs_orr[r, f], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r],
    n_fixed_bs_pfs_conf_resp_predicted[r, f], n_fixed_bs_pfs_pfs_predicted[r, f],
    summarize = FALSE
  ) |> 
    select(!c(n_bs_sample, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving))
}
  
get_sample_maturity_rvar <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      rv = map2(
        fit, stan_data, 
        \(f, d) spread_rvars(f, n_sample[l, p], maturity_rate[l, p]) |>
          bind_cols(expand.grid(d[c("lambda", "pred_week")]))
      )
    ) |> 
    unnest(rv)
}

get_all_pfs_crcr_trial_lambda_residual <- function(res) {
  spread_rvars(res, log_trial_lambda_residual[trial, t]) |> 
    mutate(trial_lambda_residual = exp(log_trial_lambda_residual)) |> 
    point_interval(log_trial_lambda_residual, trial_lambda_residual, .width = c(0.5, 0.8))  
}


get_all_pfs_crcr_trial_lambda_residual_draws <- function(res, ndraws = NULL) {
  spread_rvars(res, log_trial_lambda_residual[trial, t]) |> 
    mutate(
      log_trial_lambda_residual = thin_draws(log_trial_lambda_residual),
      trial_lambda_residual = exp(log_trial_lambda_residual)
    ) |> 
    unnest_rvars() |> 
    filter(is_null(ndraws) | (.draw <= ndraws))
}

get_pfs_pred_param <- function(res) {
  res |> 
    rowwise() |> 
    transmute(trial, rv = list(get_all_pfs_pred_param(fit))) |> 
    unnest(rv)  
}

get_cr_median_pfs_draws <- function(res, ndraws = Inf) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(spread_draws(fit, bs_cr_median_pfs[r]) |>
                  filter(.draw <= ndraws) |> # I use this to make sure all r have the same .draw 
                  left_join(as_tibble(stan_data["bootstrap_cr_maturity_rates"]) |> 
                              mutate(r = seq(n())), 
                            by = "r", relationship = "many-to-one"))
    ) |> 
    unnest(rv)
}

get_all_pfs_crcr_lambda <- function(res, stan_data = NULL) {
  spread_rvars(res, log_trial_lambda[trial, t]) |> 
    mutate(
      trial_lambda = exp(log_trial_lambda), 
      trial = if (!is_null(stan_data)) factor(trial, labels = levels(stan_data$patient_trial))
    )
}

get_all_pfs_crcr_lambda_trial_intercept <- function(res) {
  spread_rvars(res, log_lambda_gp_trial_intercept[trial]) |> 
    mutate(lambda_gp_trial_intercept = exp(log_lambda_gp_trial_intercept))
}

get_crcr_pfs_pred_param <- function(res, stan_data) {
  gather_rvars(res, covar_trial_coef[trial, m], covar_effect[trial, m], tumor_stim_pop_coef[trial, m]) |> 
    mutate(.exp_value = exp(.value)) |> 
    name_coef_indices(m, trial, stan_data)
}

get_powerscaled_variables <- function(res, metadata, stan_data) {
  metadata |> 
    rowwise() |> 
    mutate(
      ps = list(
        if (alpha == 1) res else powerscale(res, alpha = alpha, component = component, variable = c("log_crcr_trial_lambda", "log_trial_lambda"))
      )
    ) |> 
    transmute(
      alpha, component,
      baseline_hazard_rvar = list(
        gather_rvars(ps, log_crcr_trial_lambda[k, trial, t], log_trial_lambda[trial, t]) |> 
          mutate(
            .exp_value = exp(.value),
            trial = factor(trial, labels = levels(stan_data$patient_trial)),
            k = factor(k, levels = 1:2, labels = c("Non-response", "Response")) 
          ) |>  
          point_interval(.value, .exp_value, .width = c(0.5, 0.8))
      ) 
    )
}


get_covar_trial_sd <- function(res, stan_data) {
  spread_rvars(res, covar_trial_sd[m]) |> 
    name_coef_indices(m, NULL, stan_data)
}

get_joint_gng_prob <- function(mpfs_res_data, pfs6_res_data, orr_res_data, mpfs_cutoffs, pfs6_cutoffs, orr_cutoffs) {
  cutoffs <- bind_rows(mpfs = mpfs_cutoffs, pfs6 = pfs6_cutoffs, orr = orr_cutoffs, .id = "endpoint")

  bind_rows(
    mpfs = select(mpfs_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_median_pfs) |> 
      mutate(endpoint_forecast_val = weeks_to_months(endpoint_forecast_val)),
    pfs6 = select(pfs6_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_pfs6),
    orr = select(orr_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_subpop_orr),
    .id = "endpoint"
  ) |> 
    left_join(cutoffs, by = "endpoint") |> 
    pivot_wider(id_cols = c(model_type, fit_type, trial), names_from = endpoint, values_from = c(lrv, tv, endpoint_forecast_val)) |> 
    mutate(
      p_tv = Pr(endpoint_forecast_val_mpfs > tv_mpfs & endpoint_forecast_val_pfs6 > tv_pfs6 & endpoint_forecast_val_orr > tv_orr), 
      p_lrv = Pr(endpoint_forecast_val_mpfs > lrv_mpfs & endpoint_forecast_val_pfs6 > lrv_pfs6 & endpoint_forecast_val_orr > lrv_orr)
    ) 
}

# CRCR posterior functions ##############

get_conf_resp_hazard_ratios <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, patient_log_crcr_hazard_ratio[k, i]) |> 
                                  mutate(patient_crcr_hazard_ratio = exp(patient_log_crcr_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_conf_resp_hazard_ratios <- function(res, stan_data, ndraws = NULL) {
  spread_rvars(res, patient_log_crcr_hazard_ratio[k, i], ndraws = ndraws) |> 
    left_join(
      as_tibble(stan_data["patient_trial"]) |> 
        transmute(trial = patient_trial, i = seq(n())), 
      by = "i", relationship = "many-to-one"
    ) |> 
    mutate(
      patient_crcr_hazard_ratio = exp(patient_log_crcr_hazard_ratio),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response")) 
    )  
}

get_all_conf_resp_hazard_ratios_bindist <- function(res, stan_data, hb) {
  get_all_conf_resp_hazard_ratios(res, stan_data) |> 
    group_by(trial, k) |> 
    reframe(t = hb[-length(hb)], bindist = rvar_sample_hist(patient_crcr_hazard_ratio, hb))  
}

get_crcr_last_cif <- function(res, stan_data) {
  res |>
    spread_rvars(cif[i, t, k]) |> 
    # mutate(cif = thin_draws(cif)) |> 
    filter(min_rank(desc(t)) == 1) |>
    point_interval(cif, .width = c(0.5, 0.8)) |> 
    left_join(as_tibble(stan_data[c("patient_trial", "patient")]) |> mutate(i = seq(n())), by = "i", relationship = "many-to-one") |> 
    mutate(k = factor(k, levels = 1:2, labels = c("Non-response", "Response")))
}

get_all_conf_resp_cif_bindist <- function(res, stan_data, which_t, hb, ndraws = NULL) {
  spread_rvars(res, cif[i, t, k], ndraws = ndraws) |> 
    filter(t %in% which_t) |> 
    left_join(as_tibble(stan_data["patient_trial"]) |> mutate(i = seq(n())), by = "i", relationship = "many-to-one") |> 
    rename(trial = patient_trial) |> 
    group_by(trial, t, k) |> 
    reframe(p = hb[-1], cif_bindist = rvar_sample_hist(cif, hb, freq = FALSE)) |> 
    mutate(k = factor(k, levels = 1:2, labels = c("Non-response", "Response")))
}

oet_conf_resp_prob <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      prob_cause_rvars = map2(
        fit, analysis_data, 
        \(f, d) spread_rvars(f, log_prob_cause[i, k]) |>
          filter(k == 2) |>
          mutate(prob_cause = exp(log_prob_cause)) |> 
          select(!k) |> 
          left_join(transmute(d, i = seq(n()), confirmed_response), by = "i", relationship = "one-to-one")
      )
    ) |> 
    unnest(prob_cause_rvars)
}

get_conf_resp_covar_param <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, crcr_covar_effect[k, covar]), .id = "trial") 
}

get_conf_resp_tumor_param <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, crcr_tumor_stim_pop_coef[m, k]), .id = "trial") 
}

get_all_conf_resp_lambda_trial_intercept <- function(res) {
  spread_rvars(res, log_crcr_lambda_gp_trial_intercept[k, trial]) |> 
    mutate(
      crcr_lambda_gp_trial_intercept = exp(log_crcr_lambda_gp_trial_intercept),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))
    )
}

get_all_conf_resp_lambda_trial_intercept_bindist <- function(res, hb) {
  get_all_conf_resp_lambda_trial_intercept(res) |> 
    group_by(k) |> 
    reframe(p = hb[-1], bindist = rvar_sample_hist(crcr_lambda_gp_trial_intercept, hb, freq = FALSE))  
}

get_all_conf_resp_lambda <- function(res, stan_data = NULL) {
  rv <- spread_rvars(res, log_crcr_trial_lambda[k, trial, t]) |> 
    mutate(
      crcr_trial_lambda = exp(log_crcr_trial_lambda), 
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))
    )
  
  if (!is_null(stan_data)) {
    rv <- rv |> 
      mutate(trial = factor(trial, labels = levels(stan_data$patient_trial)))
  }
  
  return(rv)
}

get_all_confirmed_response <- function(res, stan_data = NULL) {
  spread_rvars(
    res, rep_confirmed_response_week[i], rep_confirmed_response_censored[i], rep_confirmed_response[i],
    forecast_confirmed_response_week[i], forecast_confirmed_response_censored[i], forecast_confirmed_response[i]
  ) |> 
    mutate(trial = stan_data$patient_trial, usubjid = stan_data$patient)  
}

get_all_confirmed_response_bindist <- function(res, stan_data, hb) {
  get_all_confirmed_response(res, stan_data) |> 
    group_by(trial) |> 
    reframe(
      t = hb[-length(hb)], 
      bindist = rvar_sample_hist(rep_confirmed_response_week, hb), 
      forecast_bindist = rvar_sample_hist(forecast_confirmed_response_week, hb) 
    )
}

get_all_rep_confirmed_trial_lambda_residual <- function(res) {
  spread_rvars(res, log_crcr_trial_lambda_residual[k, trial, t]) |> 
    mutate(crcr_trial_lambda_residual = exp(log_crcr_trial_lambda_residual)) |> 
    point_interval(log_crcr_trial_lambda_residual, crcr_trial_lambda_residual, .width = c(0.5, 0.8)) |> 
    mutate(k = factor(k, levels = 1:2, labels = c("Non-response", "Response")))
}

get_all_rep_confirmed_trial_lambda_residual_draws <- function(res, ndraws = NULL) {
  spread_rvars(res, log_crcr_trial_lambda_residual[k, trial, t]) |> 
    mutate(
      log_crcr_trial_lambda_residual = thin_draws(log_crcr_trial_lambda_residual),
      crcr_trial_lambda_residual = exp(log_crcr_trial_lambda_residual),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))
    ) |>
    unnest_rvars() |> 
    filter(is_null(ndraws) | (.draw <= ndraws))
}

prepare_cmprsk_data <- function(analysis_data) {
  analysis_data |> 
    transmute(
      confirmed_response_week, 
      confirmed_response_status = factor(confirmed_response, levels = c(FALSE, TRUE), labels = c("non-response", "response")) |>
        fct_na_value_to_level("censored") |> 
        fct_relevel("censored")
    )
}

get_obs_cif_data <- function(analysis_data) {
  analysis_data |> 
    prepare_cmprsk_data() |> 
    with(cmprsk::cuminc(confirmed_response_week, confirmed_response_status, cencode = "censored")) |> 
    map_dfr(identity, .id = "outcome") |> 
    mutate(outcome = str_remove(outcome, r"{^\d+\s+}")) |> 
    rename(estimate = est)
  
    # Not using tidycmprsk because it can't handle data that has a single observed outcome
    # cuminc(Surv(confirmed_response_week, confirmed_response_status) ~ 1, .) |> 
    # tidy()
}

get_crcr_predict_cif <- function(res, analysis_data) {
  cif <- res |>
    recover_types(analysis_data[, "trial"]) |> 
    spread_rvars(log_trial_cif[k, trial, t]) |>
    mutate(
      k = factor(k, levels = 1:2, labels = c("non-response", "response")), 
      trial_cif = exp(log_trial_cif)
    )
}

get_crcr_objective_response <- function(res, stan_data) {
  spread_rvars(res, rep_confirmed_response_forced[i], log_prob_cause[k, i]) |> 
    filter(k == 2) |>
    mutate(prob_cause = exp(log_prob_cause)) |> 
    bind_cols(stan_data[c("patient_trial", "patient", "objective_response", "confirmed_response", "confirmed_response_censored")]) |> 
    rename(trial = patient_trial)
}

get_crcr_pred_param <- function(res, stan_data) {
  gather_rvars(res, crcr_covar_trial_coef[k, trial, m], crcr_covar_effect[k, trial, m], crcr_tumor_stim_pop_coef[k, trial, m]) |> 
    mutate(.exp_value = exp(.value)) |> 
    name_coef_indices(m, trial, stan_data) |> 
    mutate(k = factor(k, levels = 1:2, labels = c("Non-response", "Response")))
}

get_crcr_covar_trial_sd <- function(res, stan_data) {
  spread_rvars(res, crcr_covar_trial_sd[m]) |> 
    name_coef_indices(m, NULL, stan_data)
}

get_orr <- function(res, analysis_data) {
  recover_types(res, analysis_data) |> 
    spread_rvars(rep_trial_orr[trial], forecast_trial_orr[trial], forecast_trial_subpop_orr[trial]) 
}

# SSM ########

get_tumor_ssls_level_param <- function(res, level = c("patient", "trial"), param, rvar_extractor = lite_spread_rvars) {
  level <- rlang::arg_match(level)
  
  # Create dynamic parameter names using the level prefix
  params <- rlang::parse_exprs(str_glue("{level}_{param}[n]"))
  
  res |> 
    rvar_extractor(!!!params) 
}

get_tumor_ssls_level_param_binned <- function(res, level, param, type, breaks = seq(-1, 1, 0.1), inv_link = exp, inv_link_breaks = exp(seq(-1, 1, 0.1))) {
  get_tumor_ssls_level_param(res, level, param, lite_gather_rvars) |>
    mutate(.rs_value = inv_link(.value), fit_type = type) |> # response scale 
    group_by(.variable, fit_type) |>
    group_modify(\(d, g) bind_rows(
      bin_point_intervals(d, .rs_value, breaks = inv_link_breaks, .width = c(0.5, 0.8)),
      bin_point_intervals(d, .value, breaks = breaks, .width = c(0.5, 0.8))
    )) |> 
    ungroup()
}

# nolint end: object_usage_linter
