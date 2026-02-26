# nolint start: object_usage_linter

get_state_patients <- function(
  analysis_data,
  sample_size = 12,
  random = TRUE,
  by = NULL,
  cond = TRUE,
  slicer = if (random) slice_sample else slice_head,
  baseline_col = mmsumdiam
) {
  analysis_data |>
    mutate(i = seq(n()), selected = {{ cond }}) |>
    select(trial, i, usubjid, visit_data, selected, pfs, right_censored) |>
    unnest(visit_data) |>
    mutate(n = seq(n())) |>
    nest(visit_data = !c(trial, i, usubjid, selected, pfs, right_censored)) |>
    filter(selected) |>
    mutate(baseline_value = map_dbl(visit_data, \(v) v |> pull({{ baseline_col }}) |> first())) |>
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

get_obs_var <- function(
  res,
  var,
  patient_states_data = NULL,
  relationship = "one-to-one",
  drop_initial = FALSE
) {
  if (drop_initial && !is_null(patient_states_data)) {
    patient_states_data <- patient_states_data |>
      group_by(i) |>
      filter(min_rank(ady) > 1) |>
      mutate(n = n - first(i)) |>
      ungroup()
  }

  rvar_data <- spread_rvars(res, {{ var }})

  if (!is_null(patient_states_data)) {
    rvar_data <- right_join(
      rvar_data,
      patient_states_data,
      by = "n",
      relationship = relationship
    )
  }

  return(rvar_data)
}

get_obs_state_var <- function(
  res,
  var,
  patient_states_data = NULL,
  drop_initial = FALSE,
  transform = identity
) {
  var_expr <- expr({{ var }}[n, p])

  get_obs_var(
    res,
    !!var_expr,
    patient_states_data,
    relationship = "many-to-one",
    drop_initial = drop_initial
  ) |>
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
  get_obs_state_var(
    res,
    patient_states_data,
    obs_patient_process_noise,
    drop_initial = TRUE
  )
}

get_obs_biomarker <- function(res, patient_states_data, var = rep_patient_log_sld[n]) {
  # Extract variable name from expression to create output column name
  var_expr <- enexpr(var)
  var_name <- if (is.call(var_expr) && var_expr[[1]] == "[") {
    as.character(var_expr[[2]])  # Extract base name from indexed expression
  } else {
    as.character(var_expr)
  }
  # Remove "log_" prefix to get biomarker name (e.g., rep_patient_log_sld -> rep_patient_sld)
  biomarker_name <- sub("_log_", "_", var_name)

  get_obs_var(res, {{ var }}, patient_states_data) |>
    mutate(
      !!biomarker_name := exp(.data[[var_name]]),
      # rh = posterior::rhat(.data[[var_name]]),
      # ess_b = posterior::ess_bulk(.data[[var_name]]),
      # ess_t = posterior::ess_tail(.data[[var_name]])
    )
}

# Backward compatibility alias
get_sld <- get_obs_biomarker

prepare_recist_data <- function(recist_rvar_data, var) {
  get_recist_simplex <- function(r, v) tibble(!!r := Pr(v == r))

  recist_rvar_data |>
    mutate(
      {{ var }} := rvar_factor(
        {{ var }},
        levels = 1:4,
        labels = c("CR", "PR", "SD", "PD")
      ),
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

get_subsample_forecast_data <- function(
  analysis_data,
  patient_states_data,
  forecast_extent = 0,
  baseline_col = mmsumdiam
) {
  max_obs_visit <- analysis_data |>
    unnest(visit_data) |>
    pull(week) |>
    max()
  overall_max_t <- max(max_obs_visit + 1, forecast_extent)

  analysis_data |>
    mutate(
      i = seq(n()),
      baseline_value = map_dbl(visit_data, \(v) v |> pull({{ baseline_col }}) |> first()),
      actual_patient_max_t = map_int(visit_data, \(d) max(d$week)),
      n_forecast_visits = overall_max_t - actual_patient_max_t
    ) |>
    filter(n_forecast_visits > 0) |>
    rowwise() |>
    reframe(
      trial,
      i,
      usubjid,
      actual_patient_max_t,
      patient_max_t,
      n_forecast_visits,
      baseline_value,
      week = seq(actual_patient_max_t + 1, overall_max_t)
    ) |>
    mutate(n = seq(n())) |>
    semi_join(patient_states_data, by = c("trial", "usubjid"))
}

get_forecast_var <- function(
  res,
  var,
  analysis_data,
  patient_states_data = analysis_data,
  relationship = "one-to-one",
  forecast_extent = 0,
  ndraws = NULL,
  baseline_col = mmsumdiam
) {
  subsample_forecast_data <- get_subsample_forecast_data(
    analysis_data,
    patient_states_data,
    forecast_extent,
    baseline_col = {{ baseline_col }}
  )

  spread_rvars(res, {{ var }}, ndraws = ndraws) |>
    right_join(subsample_forecast_data, by = "n", relationship = relationship)
}

get_forecast_state_var <- function(
  res,
  var,
  analysis_data,
  patient_states_data,
  transform = identity,
  forecast_extent = 0,
  ndraws = NULL,
  baseline_col = mmsumdiam
) {
  var_expr <- expr({{ var }}[n, p])

  get_forecast_var(
    res,
    !!var_expr,
    analysis_data,
    patient_states_data,
    relationship = "many-to-one",
    forecast_extent = forecast_extent,
    ndraws = ndraws,
    baseline_col = {{ baseline_col }}
  ) |>
    mutate(
      {{ var }} := transform({{ var }}),
      p = factor(p, levels = 1:3, labels = c("regress", "grow", "sum")),
      # rh = posterior::rhat({{ var }}),
      # ess_b = posterior::ess_bulk({{ var }}),
      # ess_t = posterior::ess_tail({{ var }})
    )
}

get_forecast_states <- function(
  res,
  analysis_data,
  patient_states_data,
  forecast_extent = 0,
  baseline_col = mmsumdiam
) {
  get_forecast_state_var(
    res,
    forecast_patient_states,
    analysis_data,
    patient_states_data,
    transform = exp,
    forecast_extent = forecast_extent,
    baseline_col = {{ baseline_col }}
  ) |>
    add_states_sum(forecast_patient_states)
}

# get_forecast_process_noise <- function(res, patient_states_data, analysis_data, forecast_extent = 0, ndraws = NULL) {
#   get_forecast_state_var(res, patient_states_data, analysis_data, forecast_patient_process_noise, forecast_extent = forecast_extent, ndraws = ndraws)
# }

get_forecast_biomarker <- function(
  res,
  analysis_data,
  patient_states_data = analysis_data,
  forecast_extent = 0,
  var = forecast_patient_log_sld[n],
  baseline_col = mmsumdiam
) {
  # Extract variable name from expression to create output column name
  var_expr <- enexpr(var)
  var_name <- if (is.call(var_expr) && var_expr[[1]] == "[") {
    as.character(var_expr[[2]])  # Extract base name from indexed expression
  } else {
    as.character(var_expr)
  }
  # Remove "log_" prefix to get biomarker name
  biomarker_name <- sub("_log_", "_", var_name)

  get_forecast_var(
    res,
    {{ var }},
    analysis_data,
    patient_states_data,
    forecast_extent = forecast_extent,
    baseline_col = {{ baseline_col }}
  ) |>
    mutate(
      !!biomarker_name := exp(.data[[var_name]]),
      # rh = posterior::rhat(.data[[var_name]]),
      # ess_b = posterior::ess_bulk(.data[[var_name]]),
      # ess_t = posterior::ess_tail(.data[[var_name]])
    )
}

# Backward compatibility alias
get_forecast_sld <- get_forecast_biomarker

get_forecast_recist <- function(
  res,
  analysis_data,
  patient_states_data = analysis_data,
  forecast_extent = 0,
  ndraws = NULL,
  baseline_col = mmsumdiam
) {
  subsample_forecast_data <- get_subsample_forecast_data(
    analysis_data,
    patient_states_data,
    forecast_extent = forecast_extent,
    baseline_col = {{ baseline_col }}
  )

  spread_rvars(res, forecast_recist[n], ndraws = ndraws) |>
    right_join(
      subsample_forecast_data,
      by = "n",
      relationship = "one-to-one"
    ) |>
    # mutate(rh = posterior::rhat(forecast_recist), ess_b = posterior::ess_bulk(forecast_recist), ess_t = posterior::ess_tail(forecast_recist)) |>
    prepare_recist_data(forecast_recist)
}

bin_point_intervals <- function(data, dist, breaks, ...) {
  data |>
    bin_dist({{ dist }}, breaks = breaks) |>
    point_interval({{ dist }}, ...)
}

# nolint end: object_usage_linter
