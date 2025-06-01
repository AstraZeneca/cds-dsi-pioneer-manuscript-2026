get_state_patients <- function(analysis_data, sample_size = 12, random = TRUE, by = NULL, cond = TRUE) {
  slicer <- if (random) slice_sample else slice
  
  analysis_data |> 
    mutate(i = seq(n()), selected = {{ cond }}) |> 
    select(trial, i, usubjid, visit_data, patient_max_t, selected) |> 
    unnest(visit_data) |> 
    mutate(n = seq(n())) |> 
    nest(visit_data = !c(trial, i, usubjid, patient_max_t, selected)) |> 
    filter(selected) |> 
    mutate(base_sld = map_dbl(visit_data, \(v) first(v$mmsumdiam))) |> 
    group_by({{ by }}) |> 
    slicer(n = sample_size) |>
    ungroup() |>
    unnest(visit_data)
}

add_states_sum <- function(states_data, states_col) {
  states_data |> 
    group_by(across(!c(p, {{ states_col }}))) |> 
    summarize(across(ends_with("states"), rvar_sum), .groups = "drop") |> 
    mutate(p = factor(3, levels = 1:3, labels = c("regress", "grow", "sum"))) |> 
    bind_rows(states_data) 
}

get_obs_var <- function(res, patient_states_data, var, drop_initial = FALSE) {
  if (drop_initial) {
      patient_states_data <- patient_states_data |> 
        group_by(i) |>
        filter(min_rank(ady) > 1) |>
        mutate(n = n - first(i)) |>
        ungroup()
  }
  
  rvar_data <- spread_rvars(res, {{ var }}) |> 
    inner_join(patient_states_data, by = "n")
}

get_obs_state_var <- function(res, patient_states_data, var, drop_initial = FALSE, transform = identity) {
  var_expr <- expr({{ var }}[n,p])
  
  get_obs_var(res, patient_states_data, !!var_expr, drop_initial) |> 
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum"))
    )
}

get_states <- function(res, patient_states_data) {
  get_obs_state_var(res, patient_states_data, states, transform = exp) |> 
    add_states_sum(states)
}

get_process_noise <- function(res, patient_states_data) {
  get_obs_state_var(res, patient_states_data, obs_patient_process_noise, drop_initial = TRUE)
}

get_sld <- function(res, patient_states_data) {
  get_obs_var(res, patient_states_data, rep_patient_log_sld[n]) |> 
    mutate(rep_patient_sld = exp(rep_patient_log_sld))
}

prepare_recist_data <- function(recist_rvar_data, var) {
  get_recist_simplex <- function(r, v) tibble(!!r := Pr(v == r))

  recist_rvar_data |>  
    mutate(
      {{ var }} := rvar_factor({{ var }}, levels = 1:4, labels = c("CR", "PR", "SD", "PD")),
      map_dfr({{ var }}, \(v) map_dfc(levels(v), \(r) get_recist_simplex(r, v))) 
    )
}


get_recist <- function(res, patient_states_data) {
  get_obs_var(res, patient_states_data, rep_recist[n]) |> 
    prepare_recist_data(rep_recist)
}

get_subsample_forecast_data <- function(patient_states_data, analysis_data, forecast_extent = 0) {
  overall_max_t <- max(max(analysis_data$patient_max_t), forecast_extent)
  
  analysis_data |> 
    mutate(
      i = seq(n()), 
      base_sld = map_dbl(visit_data, \(v) first(v$mmsumdiam)),
      n_forecast_visits = overall_max_t - patient_max_t
    ) |> 
    filter(n_forecast_visits > 0) |> 
    rowwise() |> 
    reframe(trial, i, usubjid, patient_max_t, n_forecast_visits, base_sld, week = seq(patient_max_t + 1, overall_max_t)) |> 
    mutate(n = seq(n())) |> 
    semi_join(patient_states_data, by = c("trial", "usubjid"))
}

get_forecast_var <- function(res, patient_states_data, analysis_data, var, forecast_extent = 0, ndraws = NULL) {
  subsample_forecast_data <- get_subsample_forecast_data(patient_states_data, analysis_data, forecast_extent) 
  
  spread_rvars(res, {{ var }}, ndraws = ndraws) |> 
    inner_join(subsample_forecast_data, by = "n") 
}

get_forecast_state_var <- function(res, patient_states_data, analysis_data, var, transform = identity, forecast_extent = 0, ndraws = NULL) {
  var_expr <- expr({{ var }}[n,p])
 
  get_forecast_var(res, patient_states_data, analysis_data, !!var_expr, forecast_extent = forecast_extent, ndraws = ndraws) |>  
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum"))
    )
}

get_forecast_states <- function(res, patient_states_data, analysis_data, forecast_extent = 0) {
  get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_states, transform = exp, forecast_extent = forecast_extent) |> 
    add_states_sum(forecast_patient_states)
}

get_forecast_process_noise <- function(res, patient_states_data, analysis_data, forecast_extent = 0, ndraws = NULL) {
  get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_process_noise, forecast_extent = forecast_extent, ndraws = ndraws) 
}

get_forecast_sld <- function(res, patient_states_data, analysis_data, forecast_extent = 0) {
  get_forecast_var(res, patient_states_data, analysis_data, forecast_patient_log_sld[n], forecast_extent = forecast_extent) |> 
    mutate(forecast_patient_sld = exp(forecast_patient_log_sld))
}

get_forecast_recist <- function(res, patient_states_data, analysis_data, forecast_extent = 0) {
  get_recist_simplex <- function(r, v) tibble(!!r := Pr(v == r))
  
  subsample_forecast_data <- get_subsample_forecast_data(patient_states_data, analysis_data, forecast_extent = forecast_extent) 
  
  spread_rvars(res, forecast_recist[n], ndraws = ndraws) |> 
    inner_join(subsample_forecast_data, by = "n") |> 
    # get_forecast_var(res, patient_states_data, analysis_data, forecast_recist[n]) |> 
    prepare_recist_data(forecast_recist)
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, ...)
}