library(tidyverse)
library(pointblank)

historical_patient_schema <- col_schema(
  studyid = "character",
  usubjid = "character",
  trtsdt = "date",
  trtedt = "date",
  treatment_end_week = "integer",
  treatment_end_day = "integer",
  calendar_week = "integer",
  calendar_day = "integer", 
  patient_min_t = "integer",
  patient_max_t = "integer",
  death = "logical",
  death_week = "integer",
  progression_before_death = "logical",
  right_censored = "logical",
  progress_week = "integer",
  pfs = "integer",
  interval_censored = "integer",
  age = "integer",
  ecogbl = "integer",
  baseline_albumin = "double",
  baseline_creatinine = "double",
  baseline_hemoglobin = "double",
  baseline_ldh = "double"
)

historical_patient_data <- read_csv(file.path("/mnt", "data", "PIONEER_2025_Historical_data", "cooked_patient_data.csv"))
historical_visit_data <- read_csv(file.path("/mnt", "data", "PIONEER_2025_Historical_data", "assessment_visit_data.csv"))

check_week_calc <- function(week, day) {
  week == ((day - 1) %/% 7) + 1
}

historical_action_level <- action_levels(warn_at = 1e-10, stop_at = 1)

historical_patient_agent <- create_agent(historical_patient_data, tbl_name = "patient_table", label = "Patient Level Data", actions = historical_action_level) |> 
  col_schema_match(historical_patient_schema, label = "Checking for expected column types") |> 
  rows_distinct(usubjid, label = "Unique patients") |> 
  col_vals_lte(trtsdt, vars(trtedt), label = "Treatment end date is after start date") |> 
  col_vals_gte(c(starts_with("treatment_end"), starts_with("calendar")), 1, label = "Date and index columns that should be strictly positive") |> 
  col_vals_not_null(c(starts_with("calendar"), matches("patient_(min|max)_t"), right_censored, interval_censored, pfs)) |> 
  col_vals_expr(expr = expr(check_week_calc(treatment_end_week, treatment_end_day)), label = "Check that treatment end week is correctly calculated") |> 
  col_vals_expr(expr = expr(check_week_calc(calendar_week, calendar_day)), label = "Check that calendar week is correctly calculated") |> 
  col_vals_lte(patient_min_t, vars(patient_max_t), label = "Check that min_t is lower or equal to max_t.") |> 
  col_vals_expr(expr = ~ (patient_last_visit - trtsdt) == treatment_end_day - 1, label = "End day is correctly calculated from dates") |> 
  col_vals_expr(expr = ~ patient_max_t - patient_min_t == patient_t_width - 1, label = "The t_width is correctly calculated") |> 
  col_vals_expr(expr = ~ xor(is.na(death_week), death), label = "Death required a death week") |> 
  col_vals_expr(expr = ~ xor(!progression_before_death, (!death | death_week >= progress_week)), label = "progression_before_death is correct") |> 
  col_vals_expr(expr = ~ !right_censored | interval_censored == 0, label = "Interval and right censoring as expected") |> 
  col_vals_expr(expr = ~ (right_censored & pfs == patient_max_t) | (!right_censored & pfs < progress_week), label = "PFS column is correctly calculated") |> 
  col_vals_between(age, 18, 100, label = "Check for expected ages") |>
  col_vals_expr(expr = ~ all(!is.na(studyid.y)), preconditions = \(d) left_join(d, historical_visit_data, by = "usubjid"), label = "All patients have at least one visit") |> 
  interrogate()

historical_visit_schema <- col_schema(
  studyid = "character",
  usubjid = "character",
  visitnum = "character",
  ady = "integer",
  day = "integer",
  week = "integer",
  mmsumdiam = "double",
  response = "character"
)
  
historical_visit_agent <- create_agent(historical_visit_data, tbl_name = "visit_table", label = "Visit Level Data", actions = historical_action_level) |>
  col_schema_match(historical_visit_schema, label = "Checking of expected column types") |> 
  rows_distinct(c(usubjid, visitnum, ady, day, week)) |>
  col_vals_not_null(c(ady, day, week, mmsumdiam), label = "Make sure these columns are not NA") |> 
  col_vals_not_equal(ady, 0, label = "ADY should never be zero") |> 
  col_vals_expr(expr = ~ (ady < 0 & day == ady + 1) | (ady > 0 & day == ady), label = "We're using week/day 0 as the last week/day in pre-screening") |> 
  col_vals_expr(expr = expr(check_week_calc(week, day)), label = "Check week calculation") |> 
  col_vals_gte(mmsumdiam, 0, label = "No negative SLD") |> 
  col_vals_expr(expr = ~ response %in% c("CR", "PR", "SD", "PD", "NE") | day <= 0, label = "All post-treatment visits must have a response") |> 
  col_vals_expr(expr = ~ all(!is.na(studyid.y)), preconditions = \(d) left_join(d, historical_patient_data, by = "usubjid"), label = "All visits have a corresponding patient") |> 
  interrogate()

create_multiagent(historical_patient_agent, historical_visit_agent)
