gen_patient_interval_properties <- function(fake_tumor_data, log_lambda, tumor_intercept, tumor_coef, settings) {
  pfs_model$functions$prepare_early_tumors_design_matrix(
    fake_tumor_data$tumor_size, settings$n_patient_tumors, settings$n_measures, sd(fake_tumor_data$tumor_size)
  ) |> 
    as_tibble() |> 
    set_names(c("tumor_size_1", "tumor_size_2")) |> 
    mutate(patient_id = rep(1:settings$n_patients, settings$n_patient_tumors)) |> 
    group_by(patient_id) |>
    summarize(tumor_covar = list(cbind(tumor_size_1, tumor_size_2))) |>
    rowwise() |> 
    reframe(
      patient_id, 
      progress_prob = pfs_model$functions$calculate_progress_linear_prob(log_lambda, tumor_intercept, tumor_coef, tumor_covar),
      hazard = pfs_model$functions$calculate_linear_hazard(log_lambda, tumor_intercept, tumor_coef, tumor_covar),
      survival = cumprod(progress_prob)
    )  
}

list_measures <- function(measures, n_measures) split(measures, rep(seq_along(n_measures), n_measures)) 

drop_missing_measures <- function(settings, missing_measures = NULL) {
  if (!is_null(missing_measures)) {
    settings %>%  
      list_modify(
        t_measure = list_measures(.$t_measure, .$n_measures - 1) |> 
          map(\(t) setdiff(t, missing_measures)),
        tumor_size = list_measures(.$tumor_size, rep(.$n_measures, .$n_patient_tumors)) |> 
          map(\(t) discard_at(t, missing_measures + 1)) |>  # The first one is actual for the baseline, t = 0. 
          unlist()
      ) %>%
      list_modify(
        n_measures = map_int(.$t_measure, length) + 1,
        t_measure = unlist(.$t_measure),
      )  
  } else {
    settings
  }
}

gen_fake_pfs_data <- function(patient_interval_data, settings) { 
  fake_data <- patient_interval_data |>
    nest(prob = !patient_id) |>  
    mutate(t_measure = with(settings, list_measures(t_measure, n_measures - 1))) %>% 
    mutate(
      map2(.$prob, .$t_measure, \(pd, t) pfs_model$functions$pfs_rng(pd$progress_prob, t)) |> 
        list_transpose() |> 
        set_names(c("interval_censored", "right_censored", "pfs", "actual_pfs")) |> 
        `!!!`(),
    ) %>%
    mutate(
      pfs_model$functions$identify_censoring(.$pfs, settings$n_measures, settings$t_measure) |> 
        set_names(c("stan_interval_censored", "stan_right_censored")) |> 
        `!!!`()
    )
  
  assertthat::assert_that(with(fake_data, all(stan_interval_censored == interval_censored)))
  assertthat::assert_that(with(fake_data, all(stan_right_censored == right_censored)))
  
  return(fake_data)
}

fit_sim_data <- function(settings, d, max_measures, ..., gen_pfs = TRUE, ignore_interval_censoring = FALSE) { 
  settings |> 
    list_modify(
      gen_pfs = gen_pfs,
      fit_data = TRUE,
      pfs = d$pfs, 
      ignore_interval_censoring = ignore_interval_censoring
    ) |> 
    pfs_model$sample(
      refresh = 0, 
      parallel_chains = 4, 
      init = \(chain_id) lst(tumor_stim_intercept = truncnorm::rtruncnorm(1, 0, mean = 0.01, sd = 0.05),
                             tumor_stim_coef = c(truncnorm::rtruncnorm(1, 0, mean = 0.01, sd = 0.1),
                                                 truncnorm::rtruncnorm(1, 0, mean = 0.005, sd = 0.01))),
      ...
    )
}

fit_simulations <- function(n, patient_interval_data, settings, ignore_interval_censoring = FALSE, fit_basename = NULL, tmp_dir = here("temp")) {
  fake_data_sim_with_ic <- tibble(sim_id = seq(n)) |> 
    rowwise() |> 
    mutate(sim_data = list(gen_fake_pfs_data(patient_interval_data, settings))) |> 
    ungroup() |> 
    mutate(
      sim_fit = furrr::future_map2(.progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
        sim_id, sim_data, 
        \(sid, sdata) fit_sim_data(
          settings, sdata, ignore_interval_censoring = ignore_interval_censoring, 
          output_basename = if (!is_null(fit_basename)) str_c(fit_basename, sid, sep = "_"), 
          output_dir = if (!is_null(fit_basename)) file.path(tmp_dir, "fit"))
      ), 
    ) 
}