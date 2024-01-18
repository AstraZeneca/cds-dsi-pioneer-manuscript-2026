gen_patient_interval_properties <- function(fake_tumor_data, log_lambda, tumor_intercept, tumor_coef, settings) {
  standardized <- pfs_model$functions$standardize_nonzero_tumor_sizes(fake_tumor_data$tumor_size)[[3]]
  
  t_measure_list <- with(settings, list_measures(t_measure, n_measures, n_patient_tumors))
  n_screening_t <- t_measure_list |> 
    list_flatten() |> 
    map_int(\(t) sum(t <= 0))
  
  pfs_model$functions$prepare_early_tumors_design_matrix(standardized, settings$n_patient_tumors, settings$n_measures, n_screening_t) |>  #, sd(fake_tumor_data$tumor_size)
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

list_measures <- function(measures, n_measures, n_patient_tumors) { 
  split(measures, rep(seq_along(n_measures), n_measures)) |>  
    split(rep(seq_along(n_patient_tumors), n_patient_tumors))
}

drop_missing_measures <- function(settings, measures = NULL, keep_only = FALSE) {
  if (!is_null(measures) && length(measures) > 0) {
    updated_settings <- if (is.list(measures)) {
      if (keep_only) {
        settings %>%  
          list_modify(
            t_measure = measures,
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map2(measures, \(t, m) keep_at(t, m - min(m) + 1)) |>  
              unlist()
          ) 
      } else {
        settings %>%  
          list_modify(
            t_measure = list_measures(.$t_measure, .$n_measures - 1) |> 
              map2(measures, \(t, m) setdiff(t, m)),
            tumor_size = list_measures(.$tumor_size, rep(.$n_measures, .$n_patient_tumors)) |> 
              map2(rep(measures, .$n_patient_tumors), \(t, m) discard_at(t, m + 1)) |>  
              unlist()
          ) 
      }
    } else {
      if (keep_only) {
        settings %>%  
          list_modify(
            t_measure = map(seq(sum(.$n_patient_tumors)), \(i) measures),
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) keep_at(t, measures - min(measures) + 1)) |>  
              unlist()
          ) 
      } else {
        settings %>%  
          list_modify(
            t_measure = with(., list_measures(t_measure, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) setdiff(t, measures)),
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) discard_at(t, measures - min(measures) + 1)) |>  
              unlist()
          ) 
      }
    }
    
    updated_settings %>%
      list_modify(
        n_measures = map_int(.$t_measure, length),
        t_measure = unlist(.$t_measure),
      )  
  } else {
    settings
  }
}

gen_fake_pfs_data <- function(patient_interval_data, settings) { 
  fake_data <- patient_interval_data |>
    nest(prob = !patient_id) |>  
    mutate(t_measure = with(settings, list_measures(t_measure, n_measures, n_patient_tumors))) %>% 
    mutate(
      map2(.$prob, .$t_measure, \(pd, t) pfs_model$functions$pfs_rng(pd$progress_prob, discard(first(t), \(x) x <= 0))) |> 
        list_transpose() |> 
        set_names(c("interval_censored", "right_censored", "pfs", "actual_pfs")) |> 
        `!!!`(),
    ) %>%
    mutate(
      pfs_model$functions$identify_censoring(.$pfs, rep_along(.$pfs, FALSE), settings$n_patient_tumors, settings$n_measures, settings$t_measure) |> 
        set_names(c("stan_interval_censored", "stan_right_censored")) |> 
        `!!!`()
    )
  
  assertthat::assert_that(with(fake_data, all(stan_interval_censored == interval_censored)))
  assertthat::assert_that(with(fake_data, all(stan_right_censored == right_censored)))
  
  return(fake_data)
}

create_pfs_initializer <- function(stan_data) {
  function(chain_id) { 
    lst(
      tumor_stim_coef = c(truncnorm::rtruncnorm(1, 0, mean = 0, sd = stan_data$tumor_stim_coef_sd[1]),
                          truncnorm::rtruncnorm(1, 0, mean = 0, sd = stan_data$tumor_stim_coef_sd[2])))
  }
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
      init = if (!settings$no_tumor_stim) function(chain_id) { 
        # init_vals <- lst(tumor_stim_intercept = truncnorm::rtruncnorm(1, 0, mean = 0.01, sd = 0.05))
        
        lst(
          tumor_stim_coef = c(truncnorm::rtruncnorm(1, 0, mean = 0.01, sd = 0.1),
                              truncnorm::rtruncnorm(1, 0, mean = 0.005, sd = 0.01)))
      },
      ...
    )
}

fit_simulations <- function(n, patient_interval_data, settings, ignore_interval_censoring = FALSE, fit_basename = NULL, tmp_dir = here("temp")) {
  fake_data_sim_with_ic <- tibble(sim_id = seq(n)) |> 
    rowwise() |> 
    mutate(sim_data = list(gen_fake_pfs_data(patient_interval_data, settings))) |> 
    ungroup() |> 
    transmute(
      sim_id,
      sim_data,
      sim_fit = furrr::future_map2(.progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
        sim_id, sim_data, 
        \(sid, sdata) fit_sim_data(
          settings, sdata, ignore_interval_censoring = ignore_interval_censoring, 
          output_basename = if (!is_null(fit_basename)) str_c(fit_basename, sid, sep = "_"), 
          output_dir = if (!is_null(fit_basename)) file.path(tmp_dir, "fit"))
      ), 
    ) 
}