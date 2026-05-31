# nolint start: object_usage_linter

#' Generic elicited priors for the publication covariate set
#'
#' The sclc `prepare_elicited_priors()` filters against a hardcoded
#' tribble of covariate names — most publication covariates (`male`, `ecog`,
#' `hgb`, `ldh_log`, `albumin`) are not in that tribble, so the design
#' matrix and prior vector end up with mismatched dimensions.
#'
#' This function returns a weakly-informative N(0, shrink_sd * 0.2) prior for
#' every column in the design matrix — no directional information.
prepare_publication_elicited_priors <- function(design_matrix, shrink_mean = 1.0, shrink_sd = 1.0) {
  covar_names <- colnames(design_matrix)
  if (length(covar_names) == 0) {
    return(tibble(coef_mean = numeric(0), coef_sd = numeric(0)))
  }
  tibble(
    name = factor(covar_names, levels = covar_names),
    coef_mean = 0.0,
    coef_sd = shrink_sd * 0.2
  ) |>
    arrange(name)
}

#' Build combined analysis data from the PUBLICATION target + historical CSVs.
#'
#' Replaces prepare_analysis_data() for the publication pipeline. The publication
#' CSVs have pre-computed calendar_day/calendar_week and no trtsdt, so we
#' cannot use the sclc version directly. Key differences:
#'   - potential_followup is set to patient_max_t (conservative: all censored
#'     patients are treated as admin-censored at their last visit)
#'   - original_right_censored is aliased from the CSV's right_censored column
#'     (the publication data has already applied the admin-censor adjustment)
prepare_publication_analysis_data <- function(
  target_patient_data,
  historical_patient_data,
  target_visit_data,
  historical_visit_data
) {
  bind_rows(target_patient_data, historical_patient_data) |>
    # Drop patients with no usable timing data: pfs and death_week both NA with
    # right_censored=FALSE indicates a corrupt source record (seen in sanofi_efc5505_crc).
    filter(!(!right_censored & is.na(pfs))) |>
    nest_join(
      bind_rows(target_visit_data, historical_visit_data),
      by = c("studyid", "usubjid"),
      name = "visit_data"
    ) |>
    # Drop patients with no screening visit (week <= 0): Stan requires >= 1
    # pre-screening measurement to anchor the tumor state-space model.
    filter(map_lgl(visit_data, \(v) any(v$week <= 0))) |>
    mutate(
      # Each trial×arm is a separate hierarchy group (5 groups: 1 experimental
      # + 4 historical controls). Keeps amgen control and experimental separate.
      trial = fct_drop(interaction(trial, arm, sep = "_")),
      # determine_pfs only on on-study visits (week >= 0): pre-baseline PD at
      # negative weeks produces negative det_interval_censored -> target_pfs < 0.
      map_dfr(visit_data, \(d) determine_pfs(filter(d, week >= 0), 1)),
      # Rename sex (already 0/1 numeric in publication CSVs) to `male` to match
      # the column name used in the existing publication fit and to bypass
      # prepare_covar_design_matrix's factor-relevel branch
      male = sex,
      # Median-impute missing covariates across all datasets
      across(
        c(age, ecog, hgb, ldh_log, albumin),
        \(x) replace_na(x, median(x, na.rm = TRUE))
      ),
      # Stubs required by classify_ms_patients / validate_ms_inputs:
      # original_right_censored: the CSV value is already the final adjusted flag
      original_right_censored = right_censored,
      # potential_followup: conservative assumption — patient was followed until
      # their last recorded visit (satisfies the >= patient_max_t invariant)
      potential_followup = patient_max_t,
      # Some source datasets leave progression_before_death NA when pfs == death_week
      # (same-week progression and death). Treat as FALSE (0->2 direct death) since
      # no prior documented progression event exists.
      progression_before_death = replace_na(progression_before_death, FALSE),
      # ms_prog_deterministic: TRUE if target-lesion PD was observed AND patient
      # progressed before death — Stan reads this to skip the stochastic hazard
      # contribution at T_01 for these patients (matches sclc convention)
      ms_prog_deterministic = as.integer(replace_na(
        map_lgl(visit_data, \(d) any(fct_match(d$det_response, "PD"), na.rm = TRUE)) &
          progression_before_death,
        FALSE
      )),
      # Adjust pfs upward by interval_censored for event patients (match sclc logic)
      pfs = if_else(!right_censored, pfs + interval_censored + 1L, pfs),
    ) |>
    (\(d) {
      d$ms_pattern <- classify_ms_patients(d)
      d
    })()
}

# nolint end: object_usage_linter
