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
      confirmed_response_week = if_else(
        confirmed_response,
        lag(week, default = NA),
        week
      ),
      confirmed_response_day = if_else(
        confirmed_response,
        lag(day, default = NA),
        day
      ),
      # Here lag by 2
      confirmed_response_interval_censored = confirmed_response_week -
        (if_else(
          confirmed_response,
          lag(week, n = 2L, default = 0),
          lag(week, default = 0)
        ) +
          1)
    )

  first_conf_week <- conf_resp_data |>
    drop_na(confirmed_response) |>
    filter(rank(week, ties.method = "first") == 1)

  lst(
    confirmed_response = if (nrow(first_conf_week) > 0) {
      pull(first_conf_week, confirmed_response)
    } else {
      NA
    },
    confirmed_response_censored = is.na(confirmed_response),
    confirmed_response_interval_censored = if (confirmed_response_censored) {
      0
    } else {
      pull(first_conf_week, confirmed_response_interval_censored)
    },
    confirmed_response_day = if (confirmed_response_censored) {
      max(response$day, na.rm = TRUE)
    } else {
      pull(first_conf_week, confirmed_response_day)
    },
    confirmed_response_week = if (confirmed_response_censored) {
      max(response$week, na.rm = TRUE)
    } else {
      pull(first_conf_week, confirmed_response_week)
    }
  )
}

prepare_tumor_stan_data <- function(
  analysis_data,
  covar_design_matrix,
  cond_group,
  pfs_quantiles,
  ...,
  group_col = "trial"
) {
  covar_matrix <- as.matrix(covar_design_matrix)
  n_covar <- ncol(covar_matrix)

  if (n_covar == 0) {
    covar_matrix <- array(numeric(0), dim = c(nrow(analysis_data), 0))
  }

  n_time_varying_covar <- 3L
  n_time_invariant_covar <- n_covar

  group_vec <- analysis_data[[group_col]]
  hierarchy <- create_hierarchy_structure(
    patient_assignments = list(trial = group_vec),
    n_patients = nrow(analysis_data)
  )

  lst(
    n_patients = nrow(analysis_data),
    n_trials = n_distinct(group_vec),
    patient_trial = group_vec,
    n_patient_visits = map_int(analysis_data$visit_data, nrow),
    t_patient_visits = unnest(analysis_data, visit_data) |> pull(week),
    t_patient_visits_day = unnest(analysis_data, visit_data) |> pull(ady),
    sum_tumor_size = unnest(analysis_data, visit_data) |>
      pull(mmsumdiam) |>
      divide_by(10),
    patient_t_width = analysis_data$patient_t_width,
    calendar_day = analysis_data$calendar_day,
    calendar_week = analysis_data$calendar_week,
    extend_max_all_t = 1,
    forecast_observation_interval = 6L,
    recist = unnest(analysis_data, visit_data) |>
      pull(response) |>
      factor(levels = c("CR", "PR", "SD", "PD", "NE")) |>
      fct_na_value_to_level(level = "NE"),
    target_recist = unnest(analysis_data, visit_data) |>
      pull(det_response) |>
      factor(levels = c("CR", "PR", "SD", "PD", "NE")) |>
      fct_na_value_to_level(level = "NE"),
    pfs = analysis_data$pfs,
    right_censored = analysis_data$right_censored,
    target_pfs = if_else(
      analysis_data$det_right_censored,
      analysis_data$det_pfs,
      analysis_data$det_pfs + analysis_data$det_interval_censored + 1L
    ),
    target_right_censored = analysis_data$det_right_censored,
    death_week = coalesce(analysis_data$death_week, 0L),
    n_time_varying_covar = n_time_varying_covar,
    n_time_invariant_covar = n_time_invariant_covar,
    covar_design_matrix = covar_matrix,
    n_covar = n_covar,
    n_cond_group = length(cond_group),
    cond_group_size = if (length(cond_group) > 0) map_int(cond_group, length) else integer(0),
    cond_group = if (length(cond_group) > 0) as.integer(unlist(cond_group)) else integer(0),
    pfs_quantiles = pfs_quantiles,
    n_pfs_quantiles = length(pfs_quantiles),
    n_forecast_patients = nrow(analysis_data),
    forecast_split_level = 0L,
    forecast_group = 0L,
    ms_split_level = 0L,
    n_ms_target_groups = 0L,
    ms_target_groups = integer(0),
    enable_student_t_hierarchy = 0L,
  ) |>
    list_modify(!!!hierarchy) |>
    list_modify(!!!serialize_sd_subhierarchy_modes(
      NULL, colnames(hierarchy$patient_level_groups), hierarchy$n_groups_per_level
    )) |>
    list_modify(...)
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
    pfs = {{ pfs_var }},
    death_week,
    calendar_week,
    calendar_day,
    right_censored,
    any_of("admin_right_censored_week"),
    interval_censored,
    patient = usubjid
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
  covar_formula,
  analysis_data,
  ...,
  include_covar = TRUE,
  scale_numeric = TRUE,
  handle_missing_covar = c("drop", "mean_impute")
) {
  handle_missing_covar <- arg_match(handle_missing_covar)
  covar_design_matrix <- array(NA, dim = c(nrow(analysis_data), 0))
  imputed_patients <- NULL

  if (include_covar) {
    incomplete_patients <- identify_incomplete_cases(
      analysis_data,
      covar_formula
    )

    covar_design_matrix <- withr::with_options(
      c(na.action = if (handle_missing_covar == "drop") na.omit else na.pass),
      modelr::model_matrix(analysis_data, covar_formula)
    ) |>
      select(!any_of("(Intercept)")) |>
      mutate(across(everything(), \(col) {
        scale(col, center = FALSE, scale = is.numeric(col) & scale_numeric)
      }))

    if (!is_empty(incomplete_patients)) {
      if (handle_missing_covar == "drop") {
        warning(
          length(incomplete_patients),
          " patients have incompelete cases and have been removed."
        )

        analysis_data <- slice(analysis_data, -incomplete_patients)
      } else {
        warning(
          length(incomplete_patients),
          " patients have incompelete cases. Missing covariates will be imputed."
        )
        stopifnot(scale_numeric)
        imputed_patients <- incomplete_patients
      }
    }

    covar_design_matrix <- bind_cols(
      covar_design_matrix,
      select(analysis_data, trial)
    ) |>
      group_by(trial) |>
      mutate(across(everything(), \(col) {
        scale(col, scale = FALSE) |> coalesce(0)
      })) |>
      ungroup() |>
      select(!trial)

    if (!is_empty(incomplete_patients) && handle_missing_covar != "drop") {
      warning(
        length(incomplete_patients),
        " patients have incompelete cases. Missing covariates will be imputed."
      )
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

      forecast_observation_interval = 6L,

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
      conf_resp_km = list(
        broom::tidy(survfit2(
          Surv(confirmed_response_week, 1 - confirmed_response_censored) ~ 1,
          data = as_tibble(stan_data[c(
            "confirmed_response_week",
            "confirmed_response_censored"
          )])
        )) |>
          transmute(
            t = time,
            s = estimate,
            n = n.risk,
            c = n.censor,
            e = n.event
          )
      ),

      conf_resp_km_calendar = list(
        broom::tidy(survfit2(
          Surv(
            confirmed_response_week + calendar_week - 1,
            1 - confirmed_response_censored
          ) ~ 1,
          data = as_tibble(stan_data[c(
            "confirmed_response_week",
            "calendar_week",
            "confirmed_response_censored"
          )])
        )) |>
          transmute(
            t = time,
            s = estimate,
            n = n.risk,
            c = n.censor,
            e = n.event
          )
      ),
    )
}

prepare_analysis_data <- function(trial_patient_data, trial_visit_data, admin_censor_buffer = 6L) {
  trial_dco_df <- trial_patient_data |>
    group_by(studyid) |>
    summarise(trial_dco_date = max(trtsdt + patient_max_t * 7), .groups = "drop")

  analysis_data <- nest_join(
    trial_patient_data,
    trial_visit_data,
    by = c("studyid", "usubjid"),
    name = "visit_data"
  ) |>
    select(-any_of("baseline_date")) |>
    mutate(
      across(where(is.character), as_factor),
      map_dfr(visit_data, \(d) determine_pfs(d, 1)),
      race = str_to_title(race),
      stage_num = as.integer(str_extract(stage, "\\d+")),
      across(any_of(c("first_line", "chemo_naive", "cpi_naive")), \(x) x == 1),

      target_only_right_censored = !map_lgl(visit_data, \(d) {
        any(fct_match(d$det_response, "PD"))
      }),
      non_target_pd = map_lgl(visit_data, \(d) {
        overall_pd_idx <- which(fct_match(d$response, "PD"))
        if (length(overall_pd_idx) == 0) return(FALSE)
        target_pd_idx <- which(fct_match(d$det_response, "PD"))
        if (length(target_pd_idx) == 0) return(TRUE)
        min(overall_pd_idx) < min(target_pd_idx)
      }),

      original_right_censored = right_censored,
      progression_before_death = if_else(original_right_censored, FALSE, progression_before_death),
      right_censored = right_censored |
        !((death & (death_week - patient_max_t <= admin_censor_buffer)) | progression_before_death),
      interval_censored = (1 - right_censored) * interval_censored,
      original_pfs = pfs,
      pfs = if_else(right_censored, patient_max_t, pfs),
      progression_event = !right_censored & progression_before_death,

      ms_final_state_full = case_when(
        right_censored ~ 0L,
        progression_before_death & death ~ 2L,
        progression_before_death & !death ~ 1L,
        death & !progression_before_death ~ 2L,
        TRUE ~ 1L
      ),
      ms_time_01_full = case_when(
        right_censored ~ as.integer(patient_max_t),
        progression_before_death ~ as.integer(pfs),
        TRUE ~ 0L
      ),
      ms_time_01_sld = if_else(
        !right_censored,
        as.integer(pfs),
        as.integer(patient_max_t)
      ),
      ms_censored_01_full = as.integer(right_censored | (!progression_before_death & death)),
      ms_censored_01_sld = as.integer(right_censored),

      ms_time_02 = case_when(
        death & !progression_before_death ~ as.integer(death_week),
        TRUE ~ as.integer(patient_max_t)
      ),
      ms_censored_02 = as.integer(!(death & !progression_before_death)),

      ms_max_time_02 = as.integer(pmin(ms_time_01_full, ms_time_02)),

      ms_time_12 = case_when(
        progression_before_death & death ~ pmax(1L, as.integer(death_week - pfs)),
        progression_before_death & !death ~ as.integer(patient_max_t - pfs),
        TRUE ~ 0L
      ),
      ms_censored_12 = as.integer(!(progression_before_death & death)),

      ms_prog_deterministic = as.integer(replace_na(
        map_lgl(visit_data, \(d) {
          any(fct_match(d$det_response, "PD"), na.rm = TRUE)
        }) & replace_na(progression_before_death, FALSE),
        FALSE
      )),

      trial_dco_date = trial_dco_df$trial_dco_date[
        base::match(studyid, trial_dco_df$studyid)
      ],
      potential_followup = as.numeric(trial_dco_date - trtsdt) / 7,
    ) |>
    mutate(
      pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
    )

  analysis_data$ms_pattern <- classify_ms_patients(analysis_data, admin_censor_buffer)

  analysis_data
}

prepare_covar_design_matrix <- function(analysis_data, covar_formula) {
  if (is_null(covar_formula)) {
    return(tibble(.rows = nrow(analysis_data)))
  }

  rec_data <- recipe(formula = covar_formula, data = analysis_data)
  rec_cols <- rec_data$var_info$variable

  rec_data |>
    step_mutate(
      across(any_of(c("race")), \(x) fct_relevel(x, "Other")),
      across(any_of(c("sex")), \(x) fct_relevel(x, "F")),
      across(any_of(c("smoker")), \(x) fct_relevel(x, "never")),
      across(any_of(c("pdl1_hi")), \(x) as.integer(x))
    ) |>
    (\(.) {
      r <- .
      if ("pdl1_hi" %in% rec_cols && "prev_lines" %in% rec_cols) {
        r <- r |> step_interact(~ pdl1_hi:prev_lines)
      }
      if ("pdl1_hi" %in% rec_cols || "prev_lines" %in% rec_cols) {
        r <- r |>
          step_dummy(
            all_nominal_predictors(),
            one_hot = FALSE,
            naming = \(var, lvl, ...) str_c(var, "_", lvl)
          )
      }
      r
    })() |>
    step_center(all_predictors()) |>
    step_scale(all_predictors()) |>
    prep() |>
    bake(new_data = NULL) |>
    select(!any_of("(Intercept)"))
}

discretize_pdl1 <- function(pdl1_data, cutoff = 50, handle_gte1 = "low") {
  cleaned <- str_trim(str_to_upper(pdl1_data))
  cleaned <- str_remove_all(cleaned, "%")
  cleaned <- str_replace_all(cleaned, ">=", "≥")
  cleaned <- str_replace_all(cleaned, "> =", "≥")
  cleaned <- str_replace_all(cleaned, "\\s*>\\s*", ">")

  operator <- case_when(
    str_detect(cleaned, "^≥") ~ "≥",
    str_detect(cleaned, "^>") ~ ">",
    str_detect(cleaned, "^<") ~ "<",
    str_detect(cleaned, "^≤") ~ "≤",
    TRUE ~ "="
  )

  numeric_value <- as.numeric(str_extract(cleaned, "\\d+(?:\\.\\d+)?"))

  range_pattern <- "^(\\d+)-(\\d+)$"
  range_matches <- str_match(cleaned, range_pattern)
  range_indices <- !is.na(range_matches[, 1])
  numeric_value[range_indices] <- as.numeric(range_matches[range_indices, 2])

  gte1_indices <- operator == "≥" & numeric_value == 1 & !is.na(numeric_value)

  case_when(
    gte1_indices & handle_gte1 == "low" ~ "Low",
    gte1_indices & handle_gte1 == "high" ~ "High",
    gte1_indices & handle_gte1 == "missing" ~ NA_character_,
    operator == "≥" & numeric_value >= cutoff & !is.na(numeric_value) ~ "High",
    operator == "≥" & numeric_value < cutoff & !is.na(numeric_value) ~ "Low",
    operator == ">" & numeric_value >= cutoff & !is.na(numeric_value) ~ "High",
    operator == ">" & numeric_value < cutoff & !is.na(numeric_value) ~ "Low",
    !is.na(numeric_value) & numeric_value >= cutoff ~ "High",
    !is.na(numeric_value) & numeric_value < cutoff ~ "Low",
    TRUE ~ NA_character_
  )
}

subsample_for_testing <- function(data, n_per_trial) {
  result <- data |>
    group_by(trial) |>
    group_modify(\(trial_data, trial_key) {
      stratified <- trial_data |>
        mutate(.stratum = str_c(
          ms_final_state_full, ms_prog_deterministic,
          ms_censored_12, interval_censored,
          sep = "_"
        ))

      mandatory <- stratified |>
        group_by(.stratum) |>
        slice_sample(n = 1) |>
        ungroup()

      remaining_n <- max(0L, n_per_trial - nrow(mandatory))
      if (remaining_n > 0) {
        pool <- stratified |>
          anti_join(mandatory, by = "usubjid")
        extras <- pool |>
          slice_sample(n = min(remaining_n, nrow(pool)))
        bind_rows(mandatory, extras) |> select(!.stratum)
      } else {
        mandatory |> select(!.stratum)
      }
    }) |>
    ungroup()

  cat(
    "TEST MODE: Subsampled", nrow(result), "patients from", nrow(data),
    "(", n_per_trial, "per trial )\n"
  )

  result
}

harmonize_smoker_status <- function(smoker_status) {
  smoker_status |>
    str_to_lower() |>
    fct_collapse(ever = c("ever", "former", "current"))
}

# nolint end: object_usage_linter
