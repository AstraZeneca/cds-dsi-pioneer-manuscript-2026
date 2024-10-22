#' Calculate progression-free survival from clinical data for each patient 
#'
#' @param progress_week Progress week 
#' @param death_week Death week 
#' @param right_censored Right censored 
#' @param patient_tumors of all the patient's tumors 
#'
#' @return The last observed/measured week before progression was detected
calc_pfs <- function(progress_week, right_censored, patient_tumors) {
  event_week <- progress_week
  
  # Get all the assessment weeks that happened before progression (if not censored). 
  pre_progress_weeks <- unnest(patient_tumors, tumor_history) |>
    distinct(week) |>
    filter(right_censored | week < event_week) |>
    pull(week)
  
  # There are a few patients who just have a single post treatment visit
  if (length(pre_progress_weeks) > 0) max(pre_progress_weeks) else NA_integer_
}


calc_confirmed_response <- function(response) {
  conf_resp_data <- response |> 
    mutate(
      confirmed_response = if_else(!xor(objective_response, lag(objective_response, default = NA)), objective_response, NA),
      confirmed_response_week = lag(week, default = NA),
      confirmed_response_interval_censored = confirmed_response_week - (lag(week, n = 2L, default = 0) + 1)
    ) 
  
  first_conf_week <- conf_resp_data |> 
    drop_na(confirmed_response) |> 
    filter(min_rank(week) == 1) 
  
  lst( 
    confirmed_response = if (nrow(first_conf_week) > 0) pull(first_conf_week, confirmed_response) else NA,
    confirmed_response_censored = is.na(confirmed_response),
    confirmed_response_interval_censored = if (confirmed_response_censored) 0 else pull(first_conf_week, confirmed_response_interval_censored), 
    confirmed_response_week = if (confirmed_response_censored) max(response$week, na.rm = TRUE) else pull(first_conf_week, confirmed_response_week)
  )
}

#' Convert analysis data into Stan list format data 
#'
#' @param analysis_data Analysis data frame. 
#'
#' @return Stan list data.
prepare_tumor_stan_data <- function(analysis_data) {
  lst(
    n_patients = nrow(analysis_data),
    n_trials = n_distinct(analysis_data$trial),
    tumor_location = unnest(analysis_data, patient_tumors) |> pull(tuloc) |> factor(),
    n_tumor_locations = nlevels(tumor_location), 
    patient_trial = factor(analysis_data$trial),
    n_patient_tumors = analysis_data$n_tumors,
    n_measures = analysis_data$n_measures |> unlist(),
    t_measure = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$week) |> unlist(),
    tumor_size = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$mmdiam / 10) |> unlist(),
  )
}

base_prepare_pfs_stan_data <- function(analysis_data, ..., pfs_var = pfs) {
  tumor_stan_data <- prepare_tumor_stan_data(analysis_data)
  pfs_data <- select(
      analysis_data, 
      pfs = {{ pfs_var }}, death_week, experiment_start_week, right_censored, interval_censored, patient = usubjid
    ) |> 
    mutate(
      death_week = if_else(right_censored, 0, death_week), # Death week is irrelevant if the data is censored
      patient = factor(patient) 
    ) 
  
  stan_data <- lst(
    fit_data = TRUE,
    use_tumor_model = FALSE,
    gen_pfs = TRUE,
    gen_interval_censored = FALSE,
    ignore_interval_censoring = FALSE,
    add_trial_level = FALSE,
    add_tumor_location_level = FALSE,
    fit_post_2nd_meaure_only = TRUE,
    
    fit_tumor_data = FALSE,
    gen_tumor_sizes = FALSE,
    predict_missing_sizes = FALSE, 
    multilevel_patient = FALSE,
    multilevel_tumor = FALSE,
    
    !!!tumor_stan_data,
    !!!pfs_data,
  ) |> 
    list_assign(...)
}

prepare_confirmed_resp_stan_data <- function(covar_formula, analysis_data, ..., include_covar = TRUE, scale_numeric = TRUE) {
  incomplete_patients <- c() 
  
  covar_design_matrix <- if (include_covar) {
    incomplete_patients <- which(!complete.cases(select(analysis_data, all_of(all.vars(covar_formula)))))
    
    if (!is_empty(incomplete_patients)) {
      analysis_data <- slice(analysis_data, -incomplete_patients)
      warning(length(incomplete_patients), " patients have incompelete cases and have been removed.")
    }
    
    modelr::model_matrix(analysis_data, covar_formula) |> 
      map_dfc(\(col) scale(col, scale = is.numeric(col) & scale_numeric)) |> 
      as.matrix()
  } else {
    array(NA, dim = c(nrow(analysis_data), 0))
  }
  
  pfs_stan_data <- base_prepare_pfs_stan_data(analysis_data)
  stopifnot(pfs_stan_data$n_patients == nrow(analysis_data))
  
  n_covar <- ncol(covar_design_matrix)
  
  pfs_stan_data %>% 
    list_assign(
      add_trial_level = TRUE,
      covar_design_matrix = covar_design_matrix,
      n_covar = n_covar,
      n_tumor_covar = 2,
      time_varying_conf_resp = FALSE,
      ignore_interval_censoring = FALSE,
      
      leave_out_trial = 0,
      n_bootstrap_sample = 0,
      n_bootstrap_cr_maturity_rates = 0,
      bootstrap_cr_maturity_rates = array(dim = 0),
      n_bootstrap_pfs_maturity_rates = 0,
      bootstrap_pfs_maturity_rates = array(dim = 0),
      n_fixed_bootstrap_samples = 0, 
      
      recruit_lambda = array(dim = 0), 
      recruit_phi = 0, 
      
      objective_response = analysis_data$objective_response,
      confirmed_response = coalesce(analysis_data$confirmed_response, FALSE),
      confirmed_response_censored = analysis_data$confirmed_response_censored,
      confirmed_response_interval_censored = analysis_data$confirmed_response_interval_censored,
      confirmed_response_week = analysis_data$confirmed_response_week,
    ) |>  
    list_assign(...)
}

prepare_confirmed_resp_km <- function(stan_data) {
  stan_data |> 
    rowwise() |>
    mutate(
      # conf_resp_km = list(with(
      #   stan_data, 
      #   pfs_functions$estimate_kaplan_meier(
      #     confirmed_response_week - (1 - confirmed_response_censored), # this function expects survival time not response week so we -1 for uncensored  
      #     confirmed_response_censored, 
      #     max(confirmed_response_week)
      #   )
      # )),
      
      conf_resp_km = list(broom::tidy(survfit2(
        Surv(confirmed_response_week, 1 - confirmed_response_censored) ~ 1, 
        data = as_tibble(stan_data[c("confirmed_response_week", "confirmed_response_censored")])
      )) |> select(s = estimate, n = n.risk, c = n.censor, e = n.event)),
    
      # conf_resp_km_calendar = list(with(
      #   stan_data, 
      #   pfs_functions$estimate_kaplan_meier(
      #     confirmed_response_week + experiment_start_week - 1, 
      #     confirmed_response_censored, 
      #     max(confirmed_response_week + experiment_start_week - 1)
      #   )
      # )),
      
      conf_resp_km_calendar = list(broom::tidy(survfit2(
        Surv(confirmed_response_week + experiment_start_week - 1, 1 - confirmed_response_censored) ~ 1, 
        data = as_tibble(stan_data[c("confirmed_response_week", "experiment_start_week", "confirmed_response_censored")])
      )) |> select(s = estimate, n = n.risk, c = n.censor, e = n.event)),
    )
}
