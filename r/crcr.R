get_conf_resp_hazard_ratios <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, patient_log_crcr_hazard_ratio[i, k]) |> 
                                  mutate(patient_crcr_hazard_ratio = exp(patient_log_crcr_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_conf_resp_hazard_ratios <- function(res, stan_data, ndraws = NULL) {
  spread_rvars(res, patient_log_crcr_hazard_ratio[i, k], ndraws = ndraws) |> 
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

get_crcr_predict_cif <- function(res, stan_data = NULL, week_col = rep_confirmed_response_week, reponse_col = rep_confirmed_response) {
  get_all_confirmed_response(res, stan_data) |> 
    transmute(
      trial, i,
      confirmed_response_week = {{ week_col }}, 
      confirmed_response_status = rvar_factor({{ reponse_col }}, levels = c(2, 0:1), labels = c("censored", "non-response", "response")) 
    ) |> 
    unnest_rvars() |> 
    nest(draw_data = !c(trial, .draw)) |> 
    transmute(
      trial, .draw,
      cif = map(draw_data, 
                \(d) with(d, cmprsk::cuminc(confirmed_response_week, confirmed_response_status)) |>
                  map_dfr(identity, .id = "outcome") |> 
                  mutate(outcome = str_remove(outcome, r"{^\d+\s+}")) |> 
                  rename(estimate = est)
      )
    ) |> 
    unnest(cif)
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