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
          fit, bs_cr_median_pfs[r], bs_cr_prediction_calendar_week[r], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r]
        ) |>
          mutate(n_bs_sample_cr = n_bs_sample_cr_classified + n_bs_sample_cr_unclassified) |> 
          point_interval(
            bs_cr_median_pfs, bs_cr_prediction_calendar_week, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, n_bs_sample_cr,
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
          fit, bs_pfs_median_pfs[r], bs_pfs_prediction_calendar_week[r], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r]
        ) |>
          mutate(n_bs_sample_pfs = n_bs_sample_pfs_progressed + n_bs_sample_pfs_surviving) |> 
          point_interval(
            bs_pfs_median_pfs, bs_pfs_prediction_calendar_week, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, n_bs_sample_pfs,
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
