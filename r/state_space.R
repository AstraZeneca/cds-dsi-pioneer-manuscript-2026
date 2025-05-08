get_state_patients <- function(analysis_data, random_sample = 12) {
  analysis_data |> 
    select(usubjid, visit_data) |> 
    unnest(visit_data) |> 
    mutate(n = seq(n())) |> 
    nest(visit_data = !usubjid) |> 
    mutate(i = seq(n())) |> 
    sample_n(random_sample) |> 
    unnest(visit_data)
}

get_states <- function(res, patient_states_data) {
  states_data <- spread_rvars(res, states[n, p]) |> 
    right_join(patient_states_data, by = "n") |> 
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
    inner_join(mutate(patient_states_data, n = n - 1), by = "n") 
  
  noise_data |> 
    group_by(across(!c(p, obs_patient_process_noise))) |> 
    summarize(across(ends_with("noise"), rvar_sum), .groups = "drop") |> 
    mutate(p = 3) |> 
    bind_rows(noise_data) |> 
    mutate(p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")))
}

bin_point_intervals <- function(data, dist, breaks, .width) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, .width = .width)
}