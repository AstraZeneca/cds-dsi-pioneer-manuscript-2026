<<<<<<< HEAD
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
=======
get_state_patients <- function(analysis_data, random_sample = 12) {
  analysis_data |> 
    select(usubjid, visit_data) |> 
    unnest(visit_data) |> 
    mutate(n = seq(n())) |> 
    nest(visit_data = !usubjid) |> 
    mutate(i = seq(n())) |> 
    sample_n(random_sample) |> 
    unnest(visit_data)
>>>>>>> main
}

get_states <- function(res, patient_states_data) {
  states_data <- spread_rvars(res, states[n, p]) |> 
<<<<<<< HEAD
    right_join(
      patient_states_data |> 
        group_by(i) |> 
        filter(min_rank(ady) > 1) |> 
        mutate(n = n - first(i)) |> 
        ungroup(), 
      by = "n"
    ) |> 
=======
    right_join(patient_states_data, by = "n") |> 
>>>>>>> main
    mutate(states = exp(states)) |> 
    group_by(usubjid) |> 
    mutate(states = states * first(mmsumdiam)) |> 
    ungroup()
  
  states_data |> 
    group_by(across(!c(p, states))) |> 
    summarize(across(ends_with("states"), rvar_sum), .groups = "drop") |> 
    mutate(p = 3) |> 
    bind_rows(states_data) |> 
    mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

get_process_noise <- function(res, patient_states_data) {
  noise_data <- spread_rvars(res, obs_patient_process_noise[n, p]) |> 
<<<<<<< HEAD
    inner_join(
      # mutate(patient_states_data, n = n - 1), 
      patient_states_data |> 
        group_by(i) |> 
        filter(min_rank(ady) > 1) |> 
        mutate(n = n - first(i)) |> 
        ungroup(), 
      by = "n"
    ) 
=======
    inner_join(mutate(patient_states_data, n = n - 1), by = "n") 
>>>>>>> main
  
  noise_data |> 
    group_by(across(!c(p, obs_patient_process_noise))) |> 
    summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
    mutate(p = 3) |> 
    bind_rows(noise_data) |> 
    mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

<<<<<<< HEAD
get_forecast_process_noise <- function(res, patient_states_data, analysis_data) {
  overall_max_t <- max(analysis_data$patient_max_t)
  
  subsample_forecast_data <- analysis_data |> 
    mutate(i = seq(n())) |> 
    mutate(n_forecast_visits = overall_max_t - patient_max_t) |> 
    filter(n_forecast_visits > 0) |> 
    rowwise() |> 
    reframe(i, usubjid, patient_max_t, n_forecast_visits, week = seq(patient_max_t + 1, overall_max_t)) |> 
    mutate(n = seq(n())) |> 
    semi_join(patient_states_data, by = "usubjid") 
  
  noise_data <- spread_rvars(res, forecast_patient_process_noise[n, p]) |> 
    inner_join(subsample_forecast_data, by = "n")
  
  noise_data |>
    group_by(across(!c(p, forecast_patient_process_noise))) |> 
    summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
    mutate(p = 3) |> 
    bind_rows(noise_data) |> 
    mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, ...)
=======
bin_point_intervals <- function(data, dist, breaks, .width) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, .width = .width)
>>>>>>> main
}