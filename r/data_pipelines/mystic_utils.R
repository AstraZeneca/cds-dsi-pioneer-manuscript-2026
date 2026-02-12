# =============================================================================
# FlatIron to HISTORICAL Pipeline - Utility Functions
# =============================================================================
#
# Utility functions for data transformation and validation
# Used by: flatiron_to_historical_pipeline.R
#
# =============================================================================

library(dplyr)
library(lubridate)

# =============================================================================
# DATE UTILITIES
# =============================================================================

#' Calculate study week from two dates
#' @param event_date Date of event
#' @param reference_date Reference date (usually treatment start)
#' @return Integer week number (1-indexed)
calculate_week <- function(event_date, reference_date) {
  if (is.na(event_date) | is.na(reference_date)) {
    return(NA_integer_)
  }
  as.integer((as.numeric(event_date - reference_date)) %/% 7 + 1)
}

#' Calculate study day from two dates
#' @param event_date Date of event
#' @param reference_date Reference date (usually treatment start)
#' @return Integer day number (1-indexed)
calculate_day <- function(event_date, reference_date) {
  if (is.na(event_date) | is.na(reference_date)) {
    return(NA_integer_)
  }
  as.integer(as.numeric(event_date - reference_date) + 1)
}


# =============================================================================
# PSA RESPONSE UTILITIES
# =============================================================================

#' Categorize PSA response based on percent change from baseline
#' @param pct_change Percent change from baseline (negative = decline)
#' @return Character response category
categorize_psa_response <- function(pct_change) {
  case_when(
    is.na(pct_change) ~ NA_character_,
    pct_change <= -90 ~ "PSA90",      # >=90% decline (deep response)
    pct_change <= -50 ~ "PSA50",      # >=50% decline (confirmed response)
    pct_change < 0 ~ "DECLINE",        # Any decline
    pct_change >= 25 ~ "RISE",         # >=25% rise (PSA progression)
    TRUE ~ "STABLE"
  )
}

#' Calculate PSA percent change from baseline
#' @param current_psa Current PSA value
#' @param baseline_psa Baseline PSA value
#' @return Numeric percent change
calculate_psa_pct_change <- function(current_psa, baseline_psa) {
  if (is.na(baseline_psa) | baseline_psa <= 0 | is.na(current_psa)) {
    return(NA_real_)
  }
  100 * (current_psa - baseline_psa) / baseline_psa
}


# =============================================================================
# AGE GROUP UTILITIES
# =============================================================================

#' Categorize age into groups matching clinical trial conventions
#' @param age Numeric age
#' @return Character age group
categorize_age <- function(age) {
  case_when(
    is.na(age) ~ NA_character_,
    age < 18 ~ "<18",
    age >= 18 & age < 40 ~ "18-40",
    age >= 40 & age < 65 ~ "40-65",
    age >= 65 & age < 75 ~ "65-75",
    age >= 75 ~ ">75"
  )
}


# =============================================================================
# PROGRESSION MASKING UTILITIES
# =============================================================================

#' Apply FlatIron's +14 day masking rule for progression attribution
#' 
#' FlatIron convention: Events within 14 days of treatment start are 
#' attributed to the prior treatment line, not the current one.
#' 
#' @param progression_date Date of progression
#' @param treatment_start Date of treatment start
#' @param masking_days Number of days for masking (default 14)
#' @return Logical, TRUE if progression should be attributed to this line
apply_pd_masking <- function(progression_date, treatment_start, masking_days = 14) {
  if (is.na(progression_date) | is.na(treatment_start)) {
    return(FALSE)
  }
  progression_date >= (treatment_start + masking_days)
}


# =============================================================================
# INTERVAL CENSORING UTILITIES
# =============================================================================

#' Calculate interval censored gap
#' @param progress_week Week of progression (upper bound)
#' @param last_clean_week Last assessment week without progression (lower bound)
#' @return Integer gap in weeks
calculate_interval_gap <- function(progress_week, last_clean_week) {
  if (is.na(progress_week) | is.na(last_clean_week)) {
    return(NA_integer_)
  }
  as.integer(progress_week - last_clean_week - 1)
}


# =============================================================================
# DATA VALIDATION UTILITIES
# =============================================================================

#' Validate patient data has required columns
#' @param patient_data Patient-level dataframe
#' @return Logical, TRUE if valid
validate_patient_data <- function(patient_data) {
  required_cols <- c(
    "studyid", "usubjid", "trtsdt", "trtedt",
    "age", "sex", "ecogbl", "arm",
    "death", "os_time"
  )
  
  missing_cols <- setdiff(required_cols, names(patient_data))
  
  if (length(missing_cols) > 0) {
    warning(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    return(FALSE)
  }
  
  return(TRUE)
}

#' Validate visit data has required columns
#' @param visit_data Visit-level dataframe
#' @return Logical, TRUE if valid
validate_visit_data <- function(visit_data) {
  required_cols <- c(
    "studyid", "usubjid", "visitnum", "visit_date",
    "ady", "week"
  )
  
  missing_cols <- setdiff(required_cols, names(visit_data))
  
  if (length(missing_cols) > 0) {
    warning(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    return(FALSE)
  }
  
  return(TRUE)
}


# =============================================================================
# SUMMARY UTILITIES
# =============================================================================

#' Print summary statistics for patient data
#' @param patient_data Patient-level dataframe
print_patient_summary <- function(patient_data) {
  cat("\n=== Patient Data Summary ===\n")
  cat(sprintf("Total patients: %d\n", nrow(patient_data)))
  cat(sprintf("Deaths: %d (%.1f%%)\n", 
              sum(patient_data$death, na.rm = TRUE),
              100 * mean(patient_data$death, na.rm = TRUE)))
  cat(sprintf("Progressors: %d (%.1f%%)\n",
              sum(patient_data$progressor, na.rm = TRUE),
              100 * mean(patient_data$progressor, na.rm = TRUE)))
  cat(sprintf("Median age: %.1f years\n", median(patient_data$age, na.rm = TRUE)))
  cat(sprintf("Median OS time: %.1f weeks\n", median(patient_data$os_time, na.rm = TRUE)))
  cat("\nTreatment arms:\n")
  print(table(patient_data$arm, useNA = "ifany"))
}

#' Print summary statistics for visit data
#' @param visit_data Visit-level dataframe
print_visit_summary <- function(visit_data) {
  cat("\n=== Visit Data Summary ===\n")
  cat(sprintf("Total visits: %d\n", nrow(visit_data)))
  cat(sprintf("Unique patients: %d\n", n_distinct(visit_data$usubjid)))
  cat(sprintf("Median visits per patient: %.1f\n", 
              median(table(visit_data$usubjid))))
  cat("\nResponse distribution:\n")
  print(table(visit_data$response, useNA = "ifany"))
}
