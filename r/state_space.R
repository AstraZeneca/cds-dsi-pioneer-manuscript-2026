get_state_patients <- function(analysis_data, sample_size = 12, random = TRUE, by = NULL, cond = TRUE, slicer = if (random) slice_sample else slice_head) {
  analysis_data |> 
    mutate(i = seq(n()), selected = {{ cond }}) |> 
    select(trial, i, usubjid, visit_data, patient_max_t, selected, pfs, right_censored) |> 
    unnest(visit_data) |> 
    mutate(n = seq(n())) |> 
    nest(visit_data = !c(trial, i, usubjid, patient_max_t, selected, pfs, right_censored)) |> 
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

get_obs_var <- function(res, var, patient_states_data = NULL, relationship = "one-to-one", drop_initial = FALSE) {
  if (drop_initial && !is_null(patient_states_data)) {
      patient_states_data <- patient_states_data |> 
        group_by(i) |>
        filter(min_rank(ady) > 1) |>
        mutate(n = n - first(i)) |>
        ungroup()
  }

  rvar_data <- lite_spread_rvars(res, {{ var }})
  
  if (!is_null(patient_states_data)) {
    rvar_data <- right_join(rvar_data, patient_states_data, by = "n", relationship = relationship) 
  }
  
  return(rvar_data)
}

get_obs_state_var <- function(res, var, patient_states_data = NULL, drop_initial = FALSE, transform = identity) {
  var_expr <- expr({{ var }}[n,p])
  
  get_obs_var(res, !!var_expr, patient_states_data, relationship = "many-to-one", drop_initial = drop_initial) |> 
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")),
      # rh = posterior::rhat({{ var }}), 
      # ess_b = posterior::ess_bulk({{ var }}), 
      # ess_t = posterior::ess_tail({{ var }})
    )
}

get_states <- function(res, patient_states_data) {
  get_obs_state_var(res, states, patient_states_data, transform = exp) |> 
    add_states_sum(states)
}

get_process_noise <- function(res, patient_states_data) {
  get_obs_state_var(res, patient_states_data, obs_patient_process_noise, drop_initial = TRUE)
}

get_sld <- function(res, patient_states_data) {
  get_obs_var(res, rep_patient_log_sld[n], patient_states_data) |> 
    mutate(
      rep_patient_sld = exp(rep_patient_log_sld),
      # rh = posterior::rhat(rep_patient_log_sld), 
      # ess_b = posterior::ess_bulk(rep_patient_log_sld), 
      # ess_t = posterior::ess_tail(rep_patient_log_sld)
    )
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
  get_obs_var(res, rep_recist[n], patient_states_data) |> 
    prepare_recist_data(rep_recist) |> 
    mutate(
      # rh = posterior::rhat(rep_recist), 
      # ess_b = posterior::ess_bulk(rep_recist), 
      # ess_t = posterior::ess_tail(rep_recist)
    )
}

get_subsample_forecast_data <- function(analysis_data, patient_states_data, forecast_extent = 0) {
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

get_forecast_var <- function(res, var, analysis_data, patient_states_data = analysis_data, relationship = "one-to-one", forecast_extent = 0, ndraws = NULL) {
  subsample_forecast_data <- get_subsample_forecast_data(analysis_data, patient_states_data, forecast_extent) 
  
  lite_spread_rvars(res, {{ var }}, ndraws = ndraws) |> 
    right_join(subsample_forecast_data, by = "n", relationship = relationship) 
}

get_forecast_state_var <- function(res, var, analysis_data, patient_states_data, transform = identity, forecast_extent = 0, ndraws = NULL) {
  var_expr <- expr({{ var }}[n,p])
 
  get_forecast_var(res, !!var_expr, analysis_data, patient_states_data, relationship = "many-to-one", forecast_extent = forecast_extent, ndraws = ndraws) |>  
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")),
      # rh = posterior::rhat({{ var }}), 
      # ess_b = posterior::ess_bulk({{ var }}), 
      # ess_t = posterior::ess_tail({{ var }})
    )
}

get_forecast_states <- function(res, analysis_data, patient_states_data, forecast_extent = 0) {
  get_forecast_state_var(res, forecast_patient_states, analysis_data, patient_states_data, transform = exp, forecast_extent = forecast_extent) |> 
    add_states_sum(forecast_patient_states)
}

# get_forecast_process_noise <- function(res, patient_states_data, analysis_data, forecast_extent = 0, ndraws = NULL) {
#   get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_process_noise, forecast_extent = forecast_extent, ndraws = ndraws) 
# }

get_forecast_sld <- function(res, analysis_data, patient_states_data = analysis_data, forecast_extent = 0) {
  get_forecast_var(res, forecast_patient_log_sld[n], analysis_data, patient_states_data, forecast_extent = forecast_extent) |> 
    mutate(
      forecast_patient_sld = exp(forecast_patient_log_sld),
      # rh = posterior::rhat(forecast_patient_log_sld), 
      # ess_b = posterior::ess_bulk(forecast_patient_log_sld), 
      # ess_t = posterior::ess_tail(forecast_patient_log_sld)
    )
}

get_forecast_recist <- function(res, analysis_data, patient_states_data = analysis_data, forecast_extent = 0, ndraws = NULL) {
  subsample_forecast_data <- get_subsample_forecast_data(analysis_data, patient_states_data, forecast_extent = forecast_extent) 
  
  lite_spread_rvars(res, forecast_recist[n], ndraws = ndraws) |> 
    right_join(subsample_forecast_data, by = "n", relationship = "one-to-one") |> 
    # mutate(rh = posterior::rhat(forecast_recist), ess_b = posterior::ess_bulk(forecast_recist), ess_t = posterior::ess_tail(forecast_recist)) |> 
    prepare_recist_data(forecast_recist)
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |> 
    bin_dist({{ dist }}, breaks = breaks) |> 
    point_interval({{ dist }}, ...)
}