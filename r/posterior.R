get_pfs_conf_resp_marginal_exit_prob <- function(res) {
  res |>
    ungroup() |> 
    transmute(trial, prob_rvars = map(fit, \(f) spread_rvars(f, marginal_exit_prob[i, k, t]))) |> 
    unnest(prob_rvars)
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
  
  spread_rvars(res, trial_km_est[trial, t])
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
  
   spread_rvars(res, sim_trial_median_pfs[trial]) 
}

get_pfs_conf_resp_log_hazard_ratio <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, time_invariant_log_hazard_ratio[i, k]) |> 
                                  mutate(time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_pfs_conf_resp_hazard_ratio <- function(res, stan_data) {
  spread_rvars(res, time_invariant_log_hazard_ratio[i, k]) |>
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

get_all_pfs_pred_param <- function(res) {
  gather_rvars(res, tumor_stim_pop_coef[m], covar_effect[m], conf_resp_effect) |> 
    mutate(.exp_value = exp(.value)) |> 
    mutate(
      covar = case_when(
        fct_match(.variable, "conf_resp_effect") ~ "confirmed response",
        fct_match(.variable, "tumor_stim_pop_coef") & m == 1 ~ "baseline sum of tumor sizes",
        fct_match(.variable, "tumor_stim_pop_coef") & m == 2 ~ "first post-treatment sum of tumor sizes",
        # fct_match(.variable, "covar_effect")
        m == 1 ~ "age in [18, 40)",
        m == 2 ~ "age in [40, 65)",
        m == 3 ~ "age in [65, 75)",
        m == 4 ~ "age in [75, Inf)",
        m == 5 ~ "ecog",
        m == 6 ~ "hr status: positive",
        m == 7 ~ "prior cdk46 inhibit treatment",
        m == 8 ~ "her2 status: negative",
        m == 9 ~ "her2 status: positive"
      ) |> as_factor()
    )
}

get_all_pfs_trial_pred_param <- function(res) {
  gather_rvars(res, covar_trial_coef[m, trial]) |> 
    mutate(.exp_value = exp(.value)) |> 
    mutate(
      covar = case_when(
        m == 1 ~ "baseline sum of tumor sizes",
        m == 2 ~ "first post-treatment sum of tumor sizes",
        m == 3 ~ "age in [18, 40)",
        m == 4 ~ "age in [40, 65)",
        m == 5 ~ "age in [65, 75)",
        m == 6 ~ "age in [75, Inf)",
        m == 7 ~ "ecog",
        m == 8 ~ "hr status: positive",
        m == 9 ~ "prior cdk46 inhibit treatment",
        m == 10 ~ "her2 status: negative",
        m == 11 ~ "her2 status: positive",
        m == 12 ~ "confirmed response",
      ) |> as_factor()
    )
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

get_all_pfs_crcr_lambda <- function(res) {
  spread_rvars(res, log_lambda[t]) |> 
    mutate(lambda = exp(log_lambda)) 
}

get_all_pfs_crcr_lambda_trial_intercept <- function(res) {
  spread_rvars(res, log_lambda_gp_trial_intercept[trial]) |> 
    mutate(lambda_gp_trial_intercept = exp(log_lambda_gp_trial_intercept))
}

