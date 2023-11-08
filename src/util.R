gen_patient_interval_properties <- function(fake_tumor_data, log_lambda, settings) {
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
      progress_prob = pfs_model$functions$calculate_progress_linear_prob(log_lambda, 0.2, c(0.15, 0.05), tumor_covar),
      hazard = pfs_model$functions$calculate_linear_hazard(log_lambda, 0.2, c(0.15, 0.05), tumor_covar),
      survival = cumprod(progress_prob)
    )  
}

gen_fake_pfs_data <- function(patient_interval_data, settings) { 
  patient_interval_data |>
    nest(prob = !patient_id) |> 
    mutate(
      t_measure = with(settings, split(t_measure, rep(seq(n_patients), n_measures - 1))), # -1 because t_measures doesn't include baseline measure
      pfs_res = map2(prob, t_measure, \(pd, t) pfs_model$functions$pfs_rng(pd$progress_prob, t)),
      interval_censored = map_dbl(pfs_res, \(r) r[[1]]),
      right_censored = map_dbl(pfs_res, \(r) r[[2]]),
      pfs = map_dbl(pfs_res, \(r) r[[3]]),
      actual_pfs = map_dbl(pfs_res, \(r) r[[4]]),
    )
}

fit_sim_data <- function(settings, d, max_measures, ...) { 
  settings |> 
    list_modify(
      gen_pfs = TRUE,
      fit_data = TRUE,
      pfs = d$pfs, right_censored = d$right_censored, 
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
