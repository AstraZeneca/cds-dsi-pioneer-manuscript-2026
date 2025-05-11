get_state_patients <- function(analysis_data, sample_size = 12, random = TRUE) {
  analysis_data |> 
    select(usubjid, visit_data, patient_max_t) |> 
    unnest(visit_data) |> 
    mutate(n = seq(n())) |> 
    nest(visit_data = !c(usubjid, patient_max_t)) |> 
    mutate(i = seq(n())) %>% {  
      if (random) sample_n(., sample_size) else slice(., seq(sample_size))
    } |> 
    unnest(visit_data)  
}

add_states_sum <- function(states_data, states_col) {
  states_data |> 
    group_by(across(!c(p, {{ states_col }}))) |> 
    summarize(across(ends_with("states"), rvar_sum), .groups = "drop") |> 
    mutate(p = 3) |> 
    bind_rows(states_data) |> 
    mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

get_states <- function(res, patient_states_data) {
  states_data <- spread_rvars(res, states[n, p]) |> 
    right_join(
      patient_states_data |> 
        group_by(i) |> 
        filter(min_rank(ady) > 1) |> 
        mutate(n = n - first(i)) |> 
        ungroup(), 
      by = "n"
    ) |> 
    mutate(states = exp(states)) |> 
    group_by(usubjid) |> 
    mutate(states = states * first(mmsumdiam)) |> 
    ungroup()
  
  add_states_sum(states_data, states)
  
  # states_data |> 
  #   group_by(across(!c(p, states))) |> 
  #   summarize(across(ends_with("states"), rvar_sum), .groups = "drop") |> 
  #   mutate(p = 3) |> 
  #   bind_rows(states_data) |> 
  #   mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

get_subsample_forecast_data <- function(patient_states_data, analysis_data) {
  overall_max_t <- max(analysis_data$patient_max_t)
  
  analysis_data |> 
    mutate(i = seq(n())) |> 
    mutate(n_forecast_visits = overall_max_t - patient_max_t) |> 
    filter(n_forecast_visits > 0) |> 
    rowwise() |> 
    reframe(i, usubjid, patient_max_t, n_forecast_visits, week = seq(patient_max_t + 1, overall_max_t)) |> 
    mutate(n = seq(n())) |> 
    semi_join(patient_states_data, by = "usubjid") 
}

get_forecast_states <- function(res, patient_states_data, analysis_data) {
  subsample_forecast_data <- get_subsample_forecast_data(patient_states_data, analysis_data) 
  
  forecast_states_data <- spread_rvars(res, forecast_patient_states[n, p]) |> 
    inner_join(subsample_forecast_data, by = "n") |> 
    mutate(forecast_patient_states = exp(forecast_patient_states)) 
  
  forecast_states_data |>
    add_states_sum(forecast_patient_states)
    # 
    # group_by(across(!c(p, forecast_patient_states))) |> 
    # summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
    # mutate(p = 3) |> 
    # bind_rows(forecast_states_data) |> 
    # mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

get_process_noise <- function(res, patient_states_data) {
  spread_rvars(res, obs_patient_process_noise[n, p]) |> 
    inner_join(
      patient_states_data |> 
        group_by(i) |> 
        filter(min_rank(ady) > 1) |> 
        mutate(n = n - first(i)) |> 
        ungroup(), 
      by = "n"
    ) 
  
  # noise_data |> 
  #   group_by(across(!c(p, obs_patient_process_noise))) |> 
  #   summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
  #   mutate(p = 3) |> 
  #   bind_rows(noise_data) |> 
  #   mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

get_forecast_process_noise <- function(res, patient_states_data, analysis_data) {
  subsample_forecast_data <- get_subsample_forecast_data(patient_states_data, analysis_data) 
  
  spread_rvars(res, forecast_patient_process_noise[n, p]) |> 
    inner_join(subsample_forecast_data, by = "n")
  # 
  # noise_data |>
  #   group_by(across(!c(p, forecast_patient_process_noise))) |> 
  #   summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
  #   mutate(p = 3) |> 
  #   bind_rows(noise_data) |> 
  #   mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, ...)
}