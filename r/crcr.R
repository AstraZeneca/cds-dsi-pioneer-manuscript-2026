get_conf_resp_hazard_ratios <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, patient_log_crcr_hazard_ratio[i, k]), .id = "trial") |> 
    mutate(patient_crcr_hazard_ratio = exp(patient_log_crcr_hazard_ratio)) 
}

get_conf_resp_cif <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, cif[i, t, k]), .id = "trial") 
}

get_conf_resp_prob <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      prob_cause_rvars = map2(
        fit, analysis_data, 
        \(f, d) spread_rvars(f, prob_cause[i, k]) |>
          filter(k == 2) |> 
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
    map_dfr(\(f) spread_rvars(f, crcr_covar_effect[covar, k]), .id = "trial") 
}

get_conf_resp_tumor_param <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, crcr_tumor_stim_pop_coef[m, k]), .id = "trial") 
}
