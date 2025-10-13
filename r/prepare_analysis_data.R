# nolint start: object_usage_linter

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

#' Calculate Confirmed Response from a series of response assessments
#'
#' This function processes a data frame of response assessments over time and determines
#' the confirmed response status, timing, and associated censoring information.
#'
#' @param response A data frame containing response assessment data. 
#'   Expected columns: week, day, objective_response
#'
#' @return A list containing:
#'   \item{confirmed_response}{Logical. TRUE if a confirmed response occurred, FALSE if not, NA if censored}
#'   \item{confirmed_response_censored}{Logical. TRUE if the confirmed response status is censored}
#'   \item{confirmed_response_interval_censored}{Integer. Number of weeks of interval censoring for the confirmed response}
#'   \item{confirmed_response_day}{Integer. Day of confirmed response or last assessment if censored}
#'   \item{confirmed_response_week}{Integer. Week of confirmed response or last assessment if censored}
#'
calc_confirmed_response <- function(response) {
  conf_resp_data <- response |> 
    mutate(
      confirmed_response = case_when(
        !objective_response ~ FALSE,
        objective_response & lag(objective_response, default = NA) ~ TRUE 
      ), 
      # "confirmed response" needs two consecutive CR/PR so we need to lag by 1
      confirmed_response_week = if_else(confirmed_response, lag(week, default = NA), week),
      confirmed_response_day = if_else(confirmed_response, lag(day, default = NA), day),
      confirmed_response_interval_censored = # Here lag by 2 
        confirmed_response_week - (if_else(confirmed_response, lag(week, n = 2L, default = 0), lag(week, default = 0)) + 1)
    ) 
  
  first_conf_week <- conf_resp_data |> 
    drop_na(confirmed_response) |> 
    filter(rank(week, ties.method = "first") == 1) 
  
  lst( 
    confirmed_response = if (nrow(first_conf_week) > 0) pull(first_conf_week, confirmed_response) else NA,
    confirmed_response_censored = is.na(confirmed_response),
    confirmed_response_interval_censored = if (confirmed_response_censored) 0 else pull(first_conf_week, confirmed_response_interval_censored), 
    confirmed_response_day = if (confirmed_response_censored) max(response$day, na.rm = TRUE) else pull(first_conf_week, confirmed_response_day),
    confirmed_response_week = if (confirmed_response_censored) max(response$week, na.rm = TRUE) else pull(first_conf_week, confirmed_response_week)
  )
}

#' Convert analysis data into Stan list format data 
#'
#' @param analysis_data Analysis data frame. 
#' 
#' This is used to pass to a model the tumor specific data.
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
    n_patient_visits = map_int(analysis_data$n_measures, max), 
    t_measure = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$week) |> unlist(),
    t_day_measure = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$day) |> unlist(),
    t_patient_visits = map(analysis_data$t_measure, \(t) sort(unique(unlist(t)))) |> unlist(),
    tumor_size = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$mmdiam / 10) |> unlist(),
    sum_tumor_size = map(analysis_data$tumor_sum_size, \(ts) ts$mmsumdiam / 10) |> unlist(),
    patient_t_width = analysis_data$patient_t_width,
  )
}

#' Prepare Stan data for Progression-Free Survival (PFS) analysis
#'
#' This function prepares a list of data suitable for Stan modeling of PFS,
#' combining tumor data and PFS-specific data.
#'
#' @param analysis_data A data frame containing the analysis data
#' @param ... Additional arguments to be added to or overwrite default Stan data
#' @param pfs_var The variable in analysis_data that represents PFS. Default is 'pfs'
#'
#' @return A list containing data and control parameters for Stan PFS modeling.
#'
#' @details
#' This function combines tumor data (obtained via prepare_tumor_stan_data) with
#' PFS-specific data. It sets various control parameters for the Stan model and
#' allows for these to be overwritten or supplemented via the ... argument.
base_prepare_pfs_stan_data <- function(analysis_data, ..., pfs_var = pfs) {
  tumor_stan_data <- prepare_tumor_stan_data(analysis_data)
  pfs_data <- select(
      analysis_data, 
      pfs = {{ pfs_var }}, death_week, calendar_week, calendar_day, right_censored, any_of("admin_right_censored_week"), interval_censored, patient = usubjid
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
    pfs_ignore_interval_censoring = FALSE,
    add_trial_level = FALSE,
    add_trial_level_baseline_hazard = FALSE,
    add_trial_level_prop_hazard = FALSE,
    separate_baseline_hazard = FALSE,
    separate_prop_hazard = FALSE,
    add_tumor_location_level = FALSE,
    fit_post_2nd_meaure_only = TRUE,
    pfs_only = FALSE,
    no_prop_hazard = FALSE,
    no_tumor_effects = FALSE,
    log_lik_trial = 0,
    
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

#' Identify incomplete cases in analysis data based on a covariate formula
#'
#' This function identifies rows in the analysis data that have missing values
#' for any of the variables specified in the covariate formula.
#'
#' @param analysis_data A data frame containing the analysis data
#' @param covar_formula A formula object specifying the covariates to be checked for completeness
#'
#' @return A vector of indices corresponding to rows with incomplete data
#'
#' @details
#' The function performs the following steps:
#' 1. Selects columns from the analysis data that are specified in the covariate formula
#' 2. Converts any ordered factors to unordered factors
#' 3. Identifies cases (rows) where any of the selected variables have missing values
#' 4. Returns the indices of these incomplete cases
#'
#' @examples
#' analysis_data <- data.frame(
#'   x1 = c(1, 2, NA, 4),
#'   x2 = c("A", "B", "C", NA),
#'   x3 = ordered(c("Low", "Medium", "High", "Low"))
#' )
#' covar_formula <- ~ x1 + x2 + x3
#' incomplete_cases <- identify_incomplete_cases(analysis_data, covar_formula)
#' print(incomplete_cases)  # Should return c(3, 4)
#'
identify_incomplete_cases <- function(analysis_data, covar_formula) {
  # Handle NULL formula (no covariates case)
  if (is_null(covar_formula)) {
    return(integer(0))
  }
  
  select(analysis_data, all_of(all.vars(covar_formula))) |> 
    map_if(is.ordered, \(f) factor(f, ordered = FALSE)) |> 
    complete.cases() |> 
    not() |> 
    which()
}

#' Prepare Stan data for confirmed response analysis
#'
#' This function prepares data for Stan modeling of confirmed responses in clinical trials.
#'
#' @param covar_formula Formula specifying covariates
#' @param analysis_data Data frame containing analysis data
#' @param ... Additional arguments to be passed to or to override default Stan data
#' @param include_covar Logical, whether to include covariates (default: TRUE)
#' @param scale_numeric Logical, whether to scale numeric covariates (default: TRUE)
#' @param handle_missing_covar How to handle missing covariates: "drop" or "mean_impute"
#'
#' @return A list containing prepared data for Stan modeling
#'
#' @details
#' This function processes the input data, handles covariates, and combines it with
#' PFS data to create a comprehensive dataset for Stan modeling of confirmed responses.
#' It includes options for handling missing data and scaling numeric covariates.
#'
prepare_confirmed_resp_stan_data <- function(
    covar_formula, analysis_data, ..., include_covar = TRUE, scale_numeric = TRUE, handle_missing_covar = c("drop", "mean_impute")
) {
  handle_missing_covar <- arg_match(handle_missing_covar)
  covar_design_matrix <- array(NA, dim = c(nrow(analysis_data), 0))
  imputed_patients <- NULL
  
  if (include_covar) {
    incomplete_patients <- identify_incomplete_cases(analysis_data, covar_formula) 
    
    covar_design_matrix <- withr::with_options(
      c(na.action = if (handle_missing_covar == "drop") na.omit else na.pass), 
      modelr::model_matrix(analysis_data, covar_formula)
    ) |> 
      select(!any_of("(Intercept)")) |> 
      mutate(across(everything(), \(col) scale(col, center = FALSE, scale = is.numeric(col) & scale_numeric)))
    
    if (!is_empty(incomplete_patients)) {
      if (handle_missing_covar == "drop") {
        warning(length(incomplete_patients), " patients have incompelete cases and have been removed.")
        
        analysis_data <- slice(analysis_data, -incomplete_patients)
      } else {
        warning(length(incomplete_patients), " patients have incompelete cases. Missing covariates will be imputed.")
        stopifnot(scale_numeric)
        imputed_patients <- incomplete_patients
      }
    }
    
    covar_design_matrix <- bind_cols(covar_design_matrix, select(analysis_data, trial)) |> 
      group_by(trial) |> 
      mutate(across(everything(), \(col) scale(col, scale = FALSE) |> coalesce(0))) |> 
      ungroup() |> 
      select(!trial)
    
    if (!is_empty(incomplete_patients) && handle_missing_covar != "drop") {
      warning(length(incomplete_patients), " patients have incompelete cases. Missing covariates will be imputed.")
      stopifnot(scale_numeric)
    }
  } 
  
  pfs_stan_data <- base_prepare_pfs_stan_data(analysis_data)
  stopifnot(pfs_stan_data$n_patients == nrow(analysis_data))
  stopifnot(pfs_stan_data$n_patients == nrow(covar_design_matrix))
  
  n_covar <- ncol(covar_design_matrix)
  
  pfs_stan_data |>  
    list_assign(
      add_trial_level = TRUE,
      covar_design_matrix = covar_design_matrix,
      n_covar = n_covar,
      n_tumor_covar = 2,
      time_varying_conf_resp = FALSE,
      pfs_ignore_interval_censoring = FALSE,
      crcr_ignore_interval_censoring = FALSE,
      gen_log_lik = FALSE,
      prior_sense = FALSE,
      train_beyond_cutoff = FALSE,
      
      log_lik_trial = 0,
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
      confirmed_response_day = analysis_data$confirmed_response_day,
      confirmed_response_week = analysis_data$confirmed_response_week,
      
      orr_pop = analysis_data$orr_pop,
      
      extend_max_confresp_week = 1,
      extend_max_all_t = 1,
      
      imputed_patients = imputed_patients %||% array(dim = 0),
    ) |>  
    list_assign(...)
}

#' Prepare Kaplan-Meier estimates for confirmed response
#'
#' This function calculates Kaplan-Meier estimates for confirmed response times,
#' both in study time and calendar time.
#'
#' @param stan_data A list or data frame containing the necessary data for
#'   Kaplan-Meier estimation, including 'confirmed_response_week',
#'   'confirmed_response_censored', and 'calendar_week'.
#'
#' @return A tibble with the original data and two additional list columns:
#'   \item{conf_resp_km}{A list column containing Kaplan-Meier estimates for
#'     confirmed response in study time}
#'   \item{conf_resp_km_calendar}{A list column containing Kaplan-Meier estimates for
#'     confirmed response in calendar time}
#'
#' Each Kaplan-Meier estimate list contains:
#'   \item{t}{Time points}
#'   \item{s}{Survival probability estimates}
#'   \item{n}{Number at risk}
#'   \item{c}{Number of censored observations}
#'   \item{e}{Number of events}
#'
#' @details
#' This function uses the survfit2 function to calculate Kaplan-Meier estimates.
#' It provides estimates both in study time (weeks from start of treatment) and
#' calendar time (weeks from study start).
#'
prepare_confirmed_resp_km <- function(stan_data) {
  stan_data |> 
    rowwise() |>
    mutate(
      conf_resp_km = list(broom::tidy(survfit2(
        Surv(confirmed_response_week, 1 - confirmed_response_censored) ~ 1, 
        data = as_tibble(stan_data[c("confirmed_response_week", "confirmed_response_censored")])
      )) |> transmute(t = time, s = estimate, n = n.risk, c = n.censor, e = n.event)),
      
      conf_resp_km_calendar = list(broom::tidy(survfit2(
        Surv(confirmed_response_week + calendar_week - 1, 1 - confirmed_response_censored) ~ 1, 
        data = as_tibble(stan_data[c("confirmed_response_week", "calendar_week", "confirmed_response_censored")])
      )) |> transmute(t = time, s = estimate, n = n.risk, c = n.censor, e = n.event)),
    )
}

# nolint end: object_usage_linter