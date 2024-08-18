get_conf_resp_hazard_ratios <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, patient_log_crcr_hazard_ratio[i, k]) |> 
                                  mutate(patient_crcr_hazard_ratio = exp(patient_log_crcr_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_conf_resp_cif <- function(res, which_t, hb) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(
      \(f) spread_rvars(f, cif[i, t, k]) |> 
        filter(t %in% which_t) |> 
        group_by(t, k) |> 
        reframe(p = hb[-1], cif_bindist = rvar_sample_hist(cif, hb, freq = FALSE)) |> 
        mutate(k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))),  
      .id = "trial" 
    )
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


get_rep_confirmed_response <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |>
    map_dfr(\(f) spread_rvars(f, rep_confirmed_response_week[i], rep_confirmed_response_censored[i], rep_confirmed_response[i]), .id = "trial") 
}
