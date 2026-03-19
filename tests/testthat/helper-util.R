# tests/testthat/helper-util.R
# Pure R oracles for util.stanfunctions

r_count_positive <- function(arr) sum(arr > 0L)

r_which <- function(mask, inverse = FALSE) {
  if (inverse) which(mask <= 0L) else which(mask > 0L)
}

r_rep_each <- function(to_repeat, repeats) {
  as.integer(rep(to_repeat, each = repeats))
}

r_months_to_weeks <- function(mon) mon * 365.25 / (7 * 12)

r_calendar_date_to_study_date_scalar <- function(first, cal) cal - first + 1L

r_calendar_date_to_study_date_vec <- function(first_vec, cal_scalar) {
  as.integer(cal_scalar - first_vec + 1L)
}

r_id2idx <- function(id) as.integer(id - min(id) + 1L)
