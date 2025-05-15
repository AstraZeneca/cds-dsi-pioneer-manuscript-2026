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
    mutate(p = factor(3, levels = 1:3, labels = c("regress", "grow", "sum"))) |> 
    bind_rows(states_data) 
}

get_obs_state_var <- function(res, patient_states_data, var, transform = identity) {
  var_expr <- expr({{ var }}[n,p])
  
  noise_data <- spread_rvars(res, !!var_expr) |> 
    inner_join(
      patient_states_data |> 
        group_by(i) |> 
        filter(min_rank(ady) > 1) |> 
        mutate(n = n - first(i)) |> 
        ungroup(), 
      by = "n"
    ) |>
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum"))
    )
  
  return(noise_data)
}

get_states <- function(res, patient_states_data) {
  get_obs_state_var(res, patient_states_data, states, transform = exp) |> 
    # group_by(usubjid) |> 
    # mutate(states = states * first(mmsumdiam)) |> 
    # ungroup() |> 
    add_states_sum(states)
}

get_process_noise <- function(res, patient_states_data) {
  get_obs_state_var(res, patient_states_data, obs_patient_process_noise)
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

get_forecast_state_var <- function(res, patient_states_data, analysis_data, var, transform = identity, ndraws = NULL) {
  var_expr <- expr({{ var }}[n,p])
  
  subsample_forecast_data <- get_subsample_forecast_data(patient_states_data, analysis_data) 
  
  spread_rvars(res, !!var_expr, ndraws = ndraws) |> 
    inner_join(subsample_forecast_data, by = "n") |>
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum"))
    )
}

get_forecast_states <- function(res, patient_states_data, analysis_data) {
  get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_states, transform = exp) |> 
    add_states_sum(forecast_patient_states)
}


get_forecast_process_noise <- function(res, patient_states_data, analysis_data, ndraws = NULL) {
  get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_process_noise, ndraws = ndraws) 
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, ...)
}