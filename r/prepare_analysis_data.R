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

