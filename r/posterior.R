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

get_median_pfs_conf_resp <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, sim_median_pfs))) |> 
     ungroup() |> 
     unnest(rv)
}

get_pfs_conf_resp_log_hazard_ratio <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, time_invariant_log_hazard_ratio[i, k]) |> 
                                  mutate(time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_pfs_conf_resp_bootstrap_median_pfs <- function(res) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(
        spread_rvars(fit, bs_median_pfs[r]) |>
          point_interval(bs_median_pfs, .width = c(0.5, 0.8)) |> 
          left_join(as_tibble(stan_data[c("bootstrap_cr_maturity_rates")]) |> mutate(r = seq(n())), by = "r", relationship = "many-to-one")
        # spread_rvars(fit, bootstrap_median_pfs[p], n_bootstrap_sample[p], bootstrap_maturity_rate[p]) |>
          # point_interval(bootstrap_median_pfs, n_bootstrap_sample, bootstrap_maturity_rate, .width = 0.8) |> 
          # left_join(as_tibble(stan_data[c("recruit_lambda", "prediction_week")]) |> mutate(p = seq(n())), by = "p", relationship = "one-to-one")
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
