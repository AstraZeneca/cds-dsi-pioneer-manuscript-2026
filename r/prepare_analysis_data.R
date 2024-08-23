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

prepare_analysis_data <- function(adam_adtr, adam_adtte, adam_adsl, adam_adcm, adam_adrs, pfs_functions) {
  adam_adtr |> 
    unnest(data) |> 
    filter(
      fct_match(paramcd, "LDIAM"),
      # There are a number of assessments that show up in the data but in SDTM they are marked as "NOT DONE"
      !is.na(trdy), !is.na(aval),
      tracptfl,
      aval > 0, # Drop all zero diameter assessments; these appear to be missing measurements. 
    ) |> 
    select(studyid, idap, trial,
           usubjid, trlnkid, trlnkgrp, tuloc, tulat, trtsdt, trtedt, adt, ends_with("fl"), trdy, trtp, trta, visit, visitnum, mmdiam = aval) |>
    mutate(
      week = if_else(trdy > 0, # Is post-treatment day? 
                     trdy, 
                     trdy + 1) %/% 7, # Last pre-screening week is 0. 
      treated_week = adt >= trtsdt, # Was this a post-treatment week?
      usubjid = factor(usubjid)
    ) |> 
    arrange(studyid, usubjid, trlnkid, trdy) |>
    # For each patient-tumor row, put all the tumor assessment data in a nested column
    nest(tumor_history = c(adt, trdy, week, visit, visitnum, tsplitfl, tmergefl, mmdiam, treated_week)) |>
    mutate(
      # For some weeks, there are multiple (unscheduled) visits that might happen. I'm just going to use the average for tumor sizes within those
      # weeks.
      tumor_history = map(
        tumor_history,
        \(h) nest(h, visits = c(adt, trdy, visit, visitnum, mmdiam)) |>
          mutate(
            n_week_visits = map_int(visits, nrow),
            mmdiam = map_dbl(visits, \(v) mean(v$mmdiam))
          )
      ),                                                                                          ,
      
      n_tumor_measures = map_int(tumor_history, nrow), # How many assessments/measurements per tumor

      min_t = map_int(tumor_history, \(h) min(h$week)), # First assessment week, per tumor
      max_t = map_int(tumor_history, \(h) max(h$week)), # Last assessment week, per tumor
    ) |>
    group_by(trial) |> 
    mutate(experiment_start_week = floor(time_length(trtsdt - min(trtsdt), unit = "weeks")) + 1) |> 
    ungroup() |> 
    pack(flags = ends_with("fl")) |>
    # Each row is going to be a patient and all their tumor data are put in a nested column.
    nest(patient_tumors = !c(studyid, trial, idap, usubjid, trtsdt, trtedt, trtp, trta, experiment_start_week)) |>
    left_join(# Adding a column `n_all_tumors` that counts all the tumors per patient (irrespective of being target tumors or not). 
      adam_adtr |> 
        unnest(data) |> 
        filter(fct_match(trlnkgrp, c("NON-TARGET", "TARGET", "NEW"))) |> 
        distinct(usubjid, trlnkid) |> 
        count(usubjid, name = "n_all_tumors"),
      by = "usubjid"
    ) |> 
    rowwise() |> # Per patient
    mutate(
      patient_min_t = min(patient_tumors$min_t),
      patient_max_t = max(patient_tumors$max_t),
    ) |>
    ungroup() |>
    mutate(
      patient_t_width = patient_max_t - patient_min_t + 1, # Range between first and last assessment
      n_tumors = map_int(patient_tumors, nrow), # Tumors per patient
    ) |> 
    left_join(# Get some other time-to-event variables
      adam_adtte |> 
        unnest(data) |> 
        filter(
          (fct_match(paramcd, "PFS1") & fct_match(studyid, c("DS8201-A-U301", "DS8201-A-U302", "DS8201-A-U303"))) |
            (fct_match(paramcd, "PFS") & fct_match(studyid, "DS8201-A-U201")),
          fct_match(parqual, "CENTRAL")
        ) |> 
          select(studyid, usubjid, aval, death = dthfl, dthdt, fpddt, ltmasdt, cnsr),
      by = c("studyid", "usubjid")
    ) |> 
    mutate(
      right_censored = cnsr,
      last_assessment_week = floor(time_length(ltmasdt - trtsdt, unit = "weeks") + 1), 
      death_week = floor(time_length(dthdt - trtsdt, unit = "weeks") + 1), # Week of death
      progress_week = floor(time_length(fpddt - trtsdt, unit = "weeks") + 1), # Week of disease progression
      progress_week = case_when(death & right_censored ~ pmin(death_week, last_assessment_week, na.rm = TRUE),
                                death ~ pmin(death_week, progress_week, na.rm = TRUE), # If death, progress week is death week
                                right_censored ~ pmin(last_assessment_week, progress_week, na.rm = TRUE),
                                TRUE ~ progress_week),
      aval_week = aval %/% 7,
      
      treated = # True is getting the same dosage 
        (fct_match(trial, "Breast01") & fct_match(trtp, "5.4 mg/kg")) | 
        (fct_match(trial, c("Breast02", "Breast03", "Breast04")) & fct_match(trtp, "T-DXd")),
      
      death_week = coalesce(death_week, 0L), # if NA, set to 0.
      
      n_measures = map(patient_tumors, \(tu) tu$n_tumor_measures), # A list of the number of tumor assessments.
      t_measure = map(patient_tumors, \(tu) map(tu$tumor_history, \(tr) pull(tr, week))), # A list of tumor assessment weeks
    ) |> 
    left_join(# Getting some demographic characteristics
      unnest(adam_adsl, data) |> 
        transmute(
          studyid, usubjid, country, age, sex, race, hormonr, her2ihcn, her2ish, nreg, stratar, stratav, stratr1, stratr2, stratr3, prptmyn, 
          ecogbl = factor(ecogbl, levels = 0:5, ordered = TRUE),
          region = coalesce(region1, region2),
          age_group = cut(age, c(0, 18, 40, 65, 75, Inf), right = FALSE, ordered_result = TRUE) 
        ), 
      by = c("studyid", "usubjid"),
      relationship = "one-to-one"
    ) |>
    mutate(
      hist_visceral_disease = if_else(fct_match(trial, c("Breast02", "Breast03")), fct_match(stratr3, "Y"), NA),
      her2_status = if_else(fct_match(trial, "Breast04"), 
                            factor(stratr2, levels = 1:2, labels = c("negative", "borderline")),
                            case_when(her2ihcn <= 1 ~ "negative",
                                      her2ihcn == 2 & fct_match(her2ish, c("Examined but NE", "Not Evaluable")) ~ NA,
                                      her2ihcn == 2 & !fct_match(her2ish, c("Positive", "Amplified")) ~ "borderline",
                                      TRUE ~ "positive"))
    ) |> 
    # Get prior therapy medications per patient
    nest_join(unnest(adam_adcm, data) |> select(usubjid, cmcat, response, cmdecod), by = c("usubjid"), name = "med_data") |> 
    mutate(
      med_data = map(med_data, \(cm) mutate(cm, response = str_extract(response, "CR|NA|PR|PD|SD|UNK"))),
      prior_pertuzumab_treatment = map_lgl(med_data, \(cm) with(cm, any(fct_match(cmcat, "PRIOR CANCER SYSTEMIC THERAPY") & 
                                                                          str_detect(cmdecod, regex("PERTUZUMAB", ignore_case = TRUE))))),
      prior_cdk46_inhibit_treatment = map_lgl(med_data, \(cm) with(cm, any(
        fct_match(cmcat, "PRIOR CANCER SYSTEMIC THERAPY") & str_detect(cmdecod, regex("Palbociclib|Ribociclib|Abemaciclib", ignore_case = TRUE))
      )))
    ) |> 
    # Get visit level disease response
    nest_join(
      unnest(adam_adrs, data) |> 
        filter(
          rsacptfl, 
          fct_match(parqual, "CENTRAL"), 
          fct_match(
            param, 
            c("Overall Response", 
              "Confirmed Best Overall Response",
              "Best Overall Response", 
              "Best Overall Response (03 months)", 
              "Best Overall Response (06 months)", 
              "Best Overall Response (09 months)", 
              "Best Overall Response (12 months)")
          )
        ) |>
        transmute(
          studyid, usubjid, visit, visitnum, week = ady %/% 7, param = factor(param), paramcd = factor(paramcd), 
          response = factor(avalc, levels = c("CR", "PR", "SD", "NON-CR/NON-PD", "PD", "NE", "NED"), ordered = TRUE),
        ) %>% 
        left_join(
          filter(., fct_match(param, "Overall Response")) |> 
            arrange(week) |>
            group_by(usubjid) |> 
            transmute(
              usubjid, week, visit, param,
              objective_response = case_when(response <= "PR" ~ TRUE, response == "PD" ~ FALSE, TRUE ~ NA),
            ) |> 
            ungroup(),
          by = c("usubjid", "week", "visit", "param"),
          relationship = "one-to-one"
        ),
      by = c("studyid", "usubjid"),
      name = "disease_response"
    ) |>  
    mutate(
      # The model expects "pfs" to record the number of weeks that were progress free, so I need to find the last *assessment* week before
      # disease progress was detected.
      pfs = pmap_int(lst(progress_week, right_censored, patient_tumors), calc_pfs),
    ) |> 
    filter(!is.na(pfs) & pfs >= 0) %>% # For 3 patients we can't calculate PFS because there is no visit prior to progression visit
    mutate(
      # Calculate Objective Response Rate 
      map_dfr(
        disease_response, \(r) filter(r, fct_match(param, "Overall Response")) |>
          summarize(orr_6wk = any(week <= 6 & response <= "PR"), orr_18wk = any(week <= 18 & response <= "PR"))
      ),
      
      interval_censored = with(., pfs_functions$identify_censoring(pfs, death_week, n_tumors, unlist(n_measures), unlist(t_measure))[[1]]) |> 
        pmax(0) # Issue https://github.com/azu-oncology-rd/ods-adc-early-predict-2023/issues/36
    ) 
}

prepare_raw_pfs_analysis_data <- function(analysis_data, pfs_functions) {
  analysis_data |> 
    mutate(
      patient_tumors = map(patient_tumors, \(tu) filter(tu, min_t <= 0, max_t > 0)), # Drop tumors that don't have two assessments
      
      # Recalculate these:
      n_measures = map(patient_tumors, \(tu) tu$n_tumor_measures),
      t_measure = map(patient_tumors, \(tu) map(tu$tumor_history, \(tr) pull(tr, week))),
    ) |> 
    filter(map_int(n_measures, sum) > 0) |>
    mutate(
      # After dropping tumors, these need to be recalculated 
      patient_min_t = map_int(patient_tumors, \(tu) min(tu$min_t)),
      patient_max_t = map_int(patient_tumors, \(tu) max(tu$max_t)),
      patient_t_width = patient_max_t - patient_min_t + 1,
      n_tumors = map_int(patient_tumors, nrow),
    ) %>%
    mutate(
      interval_censored = with(., pfs_functions$identify_censoring(pfs, death_week, n_tumors, unlist(n_measures), unlist(t_measure))[[1]]) |> 
        pmax(0) # Issue https://github.com/azu-oncology-rd/ods-adc-early-predict-2023/issues/36
    ) 
}

filter_pfs_analysis_data <- function(raw_pfs_data) {
  raw_pfs_data |> 
    filter(
      !is.na(pfs), !is.na(progress_week), pfs >= 0, # No missing pfs and progress must be assured to have happened after treatment
      patient_min_t <= 0, # At least one screening assessment needed,
      patient_max_t > 0 # and at least one post treatment assessment.
    )
}

prepare_orr_analysis_data <- function(analysis_data) {
  analysis_data |> 
    filter(!is.na(pfs), pfs >= 0) %>% 
    mutate(
      wk6_subsample = patient_max_t >= 6,
      wk18_subsample = patient_max_t >=18,
    )
}

calc_confirmed_response <- function(response) {
  conf_resp_data <- filter(response, param == "Overall Response") |> 
    mutate(
      confirmed_response = if_else(!xor(objective_response, lag(objective_response, default = NA)), objective_response, NA),
      confirmed_response_week = lag(week, default = NA)
    ) 
  
  first_conf_week <- conf_resp_data |> 
    drop_na(confirmed_response) |> 
    filter(min_rank(week) == 1) 
 
  lst( 
    confirmed_response = if (nrow(first_conf_week) > 0) pull(first_conf_week, confirmed_response) else NA,
    confirmed_response_censored = is.na(confirmed_response),
    confirmed_response_week = if (confirmed_response_censored) max(conf_resp_data$week) else pull(first_conf_week, confirmed_response_week)
  )
}

prepare_confirmed_resp_analysis_data <- function(analysis_data) {
  get_confirmed_response_week <- function(r, c) { 
    if (!c) (drop_na(r, objective_response) |> pull(week) |> min()) else (filter(r, fct_match(param, "Overall Response")) |> pull(week) |> max())
  }
  
  analysis_data |> 
    mutate(pfs = pmap_int(lst(progress_week = aval_week, right_censored, patient_tumors), calc_pfs)) |> 
    filter(
      !is.na(pfs), pfs >= 0,
      fct_match(hormonr, c("NEGATIVE", "POSITIVE")),
    ) |>  
    mutate(
      combined_post_treatment_t_measure = map(t_measure, \(tl) unlist(tl) |> unique() |> keep(\(t) t > 0)),
      wk6_subsample = map_lgl(combined_post_treatment_t_measure, \(t) length(t) >= 1),
      wk12_subsample = map_lgl(combined_post_treatment_t_measure, \(t) length(t) >= 2),
      wk18_subsample = map_lgl(combined_post_treatment_t_measure, \(t) length(t) >= 3),
      
      map_dfr(disease_response, calc_confirmed_response),
    )
}

get_trial_treated_obs_km <- function(sdata, pfs_functions) {
  with(sdata, {
    list(lb = pfs, ub = pfs + interval_censored) |> 
      map_dfr(function(s) {
        pfs_functions$estimate_kaplan_meier(s, right_censored, max(s)) |>
          set_names(c("s", "n", "c", "e")) |>
          as_tibble() |> 
          mutate(t = seq(0, n() - 1))
      }, .id = "btype")
    }
  )
}

get_treated_obs_km <- function(treated_trial_pfs_stan_data, pfs_functions) { 
  map_dfr(treated_trial_pfs_stan_data, \(sdata) get_trial_treated_obs_km(sdata, pfs_functions), .id = "trial")
}

#' Prepare a user-friendly data set of the first two tumor measures. 
#'
#' @param stan_data Stan list data. 
#' @param pfs_model The cmdstanr model exposing needed functions. 
#'
#' @return Data set with tumor information
get_early_tumor_pairs <- function(stan_data, pfs_functions) {
  prep_res <- with(
    stan_data, 
    pfs_functions$prepare_early_tumors_design_matrix(
      tumor_size, n_patient_tumors, n_measures, t_measure, pfs_functions$calc_n_screening_t(n_patient_tumors, n_measures, t_measure), 2
    ) 
  )
  
  prep_res[[2]] <- exec(rbind, !!!prep_res[[2]])  
  
  exec(bind_cols, !!!prep_res) |> 
    set_colnames(c("x0", "x1", "t0", "t1")) |> 
    as_tibble() |> 
    mutate(pfs = with(stan_data, rep(pfs, n_patient_tumors)), trial = with(stan_data, rep(patient_trial, n_patient_tumors)))
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
    patient_trial = analysis_data$trial,
    n_patient_tumors = analysis_data$n_tumors,
    n_measures = analysis_data$n_measures |> unlist(),
    t_measure = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$week) |> unlist(),
    tumor_size = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$mmdiam / 10) |> unlist(),
  )
}

#' Prepare analysis data for PFS model, formating in as a list for Stan. 
#'
#' @param analysis_data Pre-prepared analysis data frame.
#' @param .tumor_priors Prior parameters relevant to tumor model.
#' @param .pfs_priors Prior parameters relevant to survival model.
#' @param ... Any other variables to pass to model.
#' @param pfs_var Name of PFS variable to use from the analysis data.
#' @param orr_var Name of ORR variable to use from the analysis data.
#'
#' @return List of variables formatted for use with Stan model.
prepare_pfs_stan_data <- function(analysis_data, .tumor_priors, .pfs_priors, pfs_functions, ..., pfs_var = pfs, orr_var = orr_6wk) {
  tumor_stan_data <- prepare_tumor_stan_data(analysis_data)
  pfs_data <- select(
      analysis_data, 
      pfs = {{ pfs_var }}, orr = {{ orr_var }}, death_week, experiment_start_week, right_censored, interval_censored, patient = usubjid
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
    
    !!!.pfs_priors,
    !!!.tumor_priors,
  ) |> 
    list_assign(...)
  
  early_tumors <- get_early_tumor_pairs(stan_data, pfs_functions) |> 
    mutate(id = seq(n())) |> 
    distinct(x0, x1, .keep_all = TRUE) |> 
    pull(id)
  
  stan_data |> list_assign(grid_tumors = early_tumors, n_grid_tumors = length(early_tumors))
}

prepare_confirmed_resp_covar_formula <- function(trials = NULL) {
  covar_formula <- ~ 0 + factor(age_group, ordered = FALSE) + factor(ecogbl, ordered = FALSE) + hormonr + prior_cdk46_inhibit_treatment
  
  if (!is_null(trials) && !any(fct_match(trials, c("Breast01", "Breast04")))) {
    covar_formula <- update(covar_formula, ~ . + hist_visceral_disease)
  }
  
  if (!is_null(trials) && !all(fct_match(trials, c("Breast04")))) {
    covar_formula <- update(covar_formula, ~ . + prior_pertuzumab_treatment)
  }
  
  if (is_null(trials) || any(fct_match(trials, c("Breast04")))) {
    covar_formula <- update(covar_formula, ~ . + her2_status)
  }
  
  return(covar_formula)
}

prepare_confirmed_resp_stan_data <- function(
    covar_formula, analysis_data, .confirmed_resp_priors, .tumor_priors, .pfs_priors, pfs_functions, ..., include_covar = TRUE
) {
  pfs_stan_data <- prepare_pfs_stan_data(analysis_data, .tumor_priors, .pfs_priors, pfs_functions)
  
  stopifnot(pfs_stan_data$n_patients == nrow(analysis_data))
  
  covar_design_matrix <- if (include_covar) {
    modelr::model_matrix(analysis_data, covar_formula) |> 
      map_dfc(\(col) scale(col, scale = FALSE)) |> 
      as.matrix()
  } else {
    array(NA, dim = c(pfs_stan_data$n_patients, 0))
  }
  
  n_covar <- ncol(covar_design_matrix)
  
  pfs_stan_data %>% 
    list_assign(
      !!!.confirmed_resp_priors,
      add_trial_level = FALSE,
      use_pfs_covar = TRUE,
      covar_design_matrix = covar_design_matrix,
      n_covar = n_covar,
      n_tumor_covar = 2,
      time_varying_conf_resp = FALSE,
      ignore_interval_censoring = FALSE,
      
      prediction_week = array(dim = 0), 
      n_bootstrap_param = 0,
      n_prediction_weeks = array(dim = 0),
      recruit_lambda = array(dim = 0), 
      recruit_phi = 0, 
      n_bootstrap_samples = 0,
      
      confirmed_response = coalesce(analysis_data$confirmed_response, FALSE),
      confirmed_response_censored = analysis_data$confirmed_response_censored,
      confirmed_response_week = analysis_data$confirmed_response_week,
    ) %>% 
    list_assign(
      crcr_covar_effect_sd = rep(.$crcr_covar_effect_sd, n_covar),
      crcr_tumor_stim_pop_coef_sd = .$crcr_tumor_stim_pop_coef_sd[1:2],
      
      covar_effect_sd = rep(.$covar_effect_sd, n_covar),
      tumor_stim_pop_coef_sd = .$tumor_stim_pop_coef_sd[1:2],
    ) |> 
    list_assign(...)
}

add_bootstrap_sample <- function(stan_data, n_samples, recruit_maturity, phi, trials) {
  eligible_param <- recruit_maturity |> 
    filter(Pr(n_sample >= 10) >= 0.8) |> 
    select(trial, lambda, pred_week)
  
  if (!is_null(trials)) {
    eligible_param <- eligible_param |> filter(trial %in% trials) 
  }
  
  stan_data |> 
    list_assign(
      n_bootstrap_samples = n_samples,
      n_bootstrap_param = nrow(eligible_param),
      prediction_week = eligible_param$pred_week, 
      recruit_lambda = eligible_param$lambda,
      recruit_phi = phi 
    )
}

get_ic_data <- function(pfs_analysis_data, tumor_priors, pfs_priors, pfs_functions) {
  pfs_analysis_data %>% 
    mutate(
      interval_censoring = prepare_pfs_stan_data(., tumor_priors, pfs_priors, pfs_functions) |>
        with(pfs_functions$identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure)) |> 
        pluck(1)
    )
}

filter_conf_resp_missing_covar <- function(analysis_data, covar_formula) {
  filter(analysis_data, if_all(all_of(all.vars(covar_formula)), \(coef) !is.na(coef)))
}

prepare_trial_confirmed_resp_stan_data <- function(analysis_data, confirmed_resp_priors, tumor_priors, pfs_priors, pfs_functions) {
  analysis_data |> 
    filter(wk12_subsample) |> 
    nest(.by = trial, .key = "analysis_data")  |> 
    rowwise() |> 
    mutate(
      covar_formula = list(prepare_confirmed_resp_covar_formula(trial)),
      # Drop rows that have missing covars 
      analysis_data = list(filter_conf_resp_missing_covar(analysis_data, covar_formula)),
      stan_data = list(
        prepare_confirmed_resp_stan_data(
          covar_formula, analysis_data, 
          confirmed_resp_priors, .tumor_priors = tumor_priors, .pfs_priors = pfs_priors, 
          pfs_functions = pfs_functions
        )
      ),
      # init_fun = list(create_pfs_initializer(stan_data)), 
      init_fun = list(if (fct_match(trial, "Breast02")) 0), 
    )
}

prepare_confirmed_resp_obs_km <- function(all_stan_data, pfs_functions) {
  all_stan_data |> 
    rowwise() |> 
    mutate(obs_km = list(get_trial_treated_obs_km(stan_data, pfs_functions)))
}

prepare_confirmed_resp_km <- function(stan_data) {
  stan_data |> 
    rowwise() |>
    mutate(
      conf_resp_km = list(with(
        stan_data, 
        pfs_functions$estimate_kaplan_meier(confirmed_response_week, confirmed_response_censored, max(confirmed_response_week))
      )),
    
      conf_resp_km_calendar = list(with(
        stan_data, 
        pfs_functions$estimate_kaplan_meier(
          confirmed_response_week + experiment_start_week - 1, 
          confirmed_response_censored, 
          max(confirmed_response_week + experiment_start_week - 1)
        )
      ))
    )
}