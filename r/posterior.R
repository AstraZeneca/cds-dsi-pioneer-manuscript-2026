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

get_all_pfs_conf_resp_log_hazard_ratio <- function(res, stan_data) {
  spread_rvars(res, time_invariant_log_hazard_ratio[i, k]) |>
    mutate(time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio)) |> 
    left_join(as_tibble(stan_data["patient_trial"]) |> mutate(i = seq(n())), by = "i", relationship = "many-to-one") |> 
    rename(trial = patient_trial)  
}

get_pfs_conf_resp_bootstrap_cr_median_pfs <- function(res) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(
        spread_rvars(
          fit, 
          bs_cr_prediction_calendar_week[r], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r],
          bs_cr_median_pfs[r], bs_cr_orr[r],
          n_bs_cr_conf_resp_predicted[r], n_bs_cr_pfs_predicted[r]
        ) |>
          mutate(
            n_bs_sample_cr = n_bs_sample_cr_classified + n_bs_sample_cr_unclassified,
            across(ends_with("predicted"), \(n)  n / n_bs_sample_cr, .names = "{.col}_prop")
          ) |> 
          unnest_rvars() |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
          point_interval(
            bs_cr_prediction_calendar_week, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, n_bs_sample_cr,
            bs_cr_median_pfs, bs_cr_orr, 
            n_bs_cr_conf_resp_predicted, n_bs_cr_pfs_predicted, n_bs_cr_conf_resp_predicted_prop, n_bs_cr_pfs_predicted_prop,
            na.rm = TRUE, # if n_bs_sample_cr is 0, we'll get some NaNs. These are few so we'll bite the bullet and drop them.
            .width = c(0.5, 0.8)
          ) |> 
          left_join(as_tibble(stan_data[c("bootstrap_cr_maturity_rates")]) |> mutate(r = seq(n())), by = "r", relationship = "many-to-one")
      )
    ) |> 
    ungroup() |> 
    unnest(rv)  
}    

get_pfs_conf_resp_bootstrap_pfs_median_pfs <- function(res) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(
        spread_rvars(
          fit, 
          bs_pfs_prediction_calendar_week[r], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r],
          bs_pfs_median_pfs[r], bs_pfs_orr[r],
          n_bs_pfs_conf_resp_predicted[r], n_bs_pfs_pfs_predicted[r]
        ) |>
          mutate(
            n_bs_sample_pfs = n_bs_sample_pfs_progressed + n_bs_sample_pfs_surviving,
            across(ends_with("predicted"), \(n) n / n_bs_sample_pfs, .names = "{.col}_prop")
          ) |> 
          unnest_rvars() |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
          point_interval(
            bs_pfs_prediction_calendar_week, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, n_bs_sample_pfs,
            bs_pfs_median_pfs, bs_pfs_orr, 
            n_bs_pfs_conf_resp_predicted, n_bs_pfs_pfs_predicted, n_bs_pfs_conf_resp_predicted_prop, n_bs_pfs_pfs_predicted_prop,
            na.rm = TRUE, # if n_bs_sample_cr is 0, we'll get some NaNs. These are few so we'll bite the bullet and drop them.
            .width = c(0.5, 0.8)
          ) |> 
          left_join(as_tibble(stan_data[c("bootstrap_pfs_maturity_rates")]) |> mutate(r = seq(n())), by = "r", relationship = "many-to-one")
      )
    ) |> 
    ungroup() |> 
    unnest(rv)  
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
      )
    )
}

get_pfs_pred_param <- function(res) {
  res |> 
    rowwise() |> 
    transmute(trial, rv = list(get_all_pfs_pred_param(fit))) |> 
    unnest(rv)  
}