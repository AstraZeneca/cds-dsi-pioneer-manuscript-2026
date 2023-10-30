gen_patient_interval_properties <- function(fake_tumor_data, log_lambda) {
  fake_tumor_data |>
    mutate(across(c(tumor_size_1, tumor_size_2), \(tsize) tsize / sd(tsize))) |> 
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

gen_fake_pfs_data <- function(patient_interval_data) { 
  patient_interval_data |> 
    group_by(patient_id) |> 
    summarize(pfs = pfs_model$functions$pfs_rng(progress_prob)) |> 
    mutate(right_censored = pfs >= length(log_lambda), interval_censored = 0)
}

fit_sim_data <- function(settings, d, ...) { 
  settings |> 
    list_modify(
      gen_pfs = TRUE,
      fit_data = TRUE,
      pfs = d$pfs, right_censored = d$right_censored
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
