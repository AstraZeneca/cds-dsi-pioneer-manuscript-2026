# tests/testthat/helper-lfo.R
#
# Pure R mirror of get_testing_visit_week_bounds() from stan/lfo.stanfunctions.
# Used to compute reference expected values for Stan tests.

#
# Pre-conditions (same as Stan asserts):
#   - t_patient_visits_week and t_patient_visits_day are strictly ascending
#     within each patient's visit block
#   - Every patient has at least 1 visit (patient_visit_pos[i+1] > patient_visit_pos[i])
#   - patient_visit_pos has length n_patients + 1


#' R mirror of Stan's get_testing_visit_week_bounds().
#'
#' @param oos_patient_idx     integer vector [n_cutoffs]
#' @param last_visit_calendar_day_sort_idx integer vector [n_patients]
#' @param cutoff_calendar_day integer vector [n_cutoffs]
#' @param patient_calendar_day integer vector [n_patients]
#' @param t_patient_visits_week integer vector [n_visits]
#' @param t_patient_visits_day  integer vector [n_visits]
#' @param patient_visit_pos     integer vector [n_patients + 1]
#'
#' @return named list:
#'   first_testing_visit_week  matrix[n_cutoffs, n_patients]
#'   last_testing_visit_week   array [n_cutoffs, n_cutoffs, n_patients]
#'   testing_start_idx         matrix[n_cutoffs, n_patients]
#'   testing_end_idx           array [n_cutoffs, n_cutoffs, n_patients]
r_get_testing_visit_week_bounds <- function(
  oos_patient_idx,
  last_visit_calendar_day_sort_idx,
  cutoff_calendar_day,
  patient_calendar_day,
  t_patient_visits_week,
  t_patient_visits_day,
  patient_visit_pos
) {
  n_patients <- length(patient_calendar_day)
  n_futures  <- length(oos_patient_idx)
  min_all_t  <- min(t_patient_visits_week)

  # Validate pre-conditions (mirrors Stan's assert_strict_ascending calls)
  for (i in seq_len(n_patients)) {
    vs <- patient_visit_pos[i]
    ve <- patient_visit_pos[i + 1L] - 1L
    stopifnot(ve >= vs)  # at least one visit
    if (ve > vs) {
      weeks <- t_patient_visits_week[vs:ve]
      days  <- t_patient_visits_day[vs:ve]
      stopifnot(all(diff(weeks) > 0), all(diff(days) > 0))
    }
  }

  first_testing_visit_week <- matrix(0L, nrow = n_futures, ncol = n_patients)
  last_testing_visit_week  <- array(min_all_t, dim = c(n_futures, n_futures, n_patients))
  testing_start_idx        <- matrix(0L, nrow = n_futures, ncol = n_patients)
  testing_end_idx          <- array(0L,  dim = c(n_futures, n_futures, n_patients))

  for (n in seq_len(n_futures)) {
    n_curr_patients <- n_patients - oos_patient_idx[n] + 1L
    curr_patients   <- last_visit_calendar_day_sort_idx[oos_patient_idx[n]:n_patients]

    # Lower cutoff: patient-specific study day for cutoff n (+1 is intentional)
    lower_sd <- cutoff_calendar_day[n] - patient_calendar_day[curr_patients] + 1L

    for (i_idx in seq_len(n_curr_patients)) {
      i           <- curr_patients[i_idx]
      visit_start <- patient_visit_pos[i]
      visit_end   <- patient_visit_pos[i + 1L] - 1L
      n_visits    <- visit_end - visit_start + 1L
      p_days      <- t_patient_visits_day[visit_start:visit_end]
      p_weeks     <- t_patient_visits_week[visit_start:visit_end]

      # First visit STRICTLY AFTER lower cutoff (mirrors Stan while loop)
      t_idx <- 1L
      while (t_idx <= n_visits && p_days[t_idx] <= max(0L, lower_sd[i_idx])) {
        t_idx <- t_idx + 1L
      }
      if (t_idx <= n_visits) {
        first_testing_visit_week[n, i] <- p_weeks[t_idx]
        testing_start_idx[n, i]        <- visit_start + t_idx - 1L
      }
    }

    # Upper cutoff: last visit on or before cutoff m — only computed for m > n
    # NOTE: guard required because R's (n+1):n gives c(n+1, n), not an empty range
    if (n < n_futures) {
      for (m in (n + 1L):n_futures) {
        upper_sd <- cutoff_calendar_day[m] - patient_calendar_day[curr_patients] + 1L

        for (i_idx in seq_len(n_curr_patients)) {
          i           <- curr_patients[i_idx]
          visit_start <- patient_visit_pos[i]
          visit_end   <- patient_visit_pos[i + 1L] - 1L
          n_visits    <- visit_end - visit_start + 1L
          rev_days    <- rev(t_patient_visits_day[visit_start:visit_end])
          rev_weeks   <- rev(t_patient_visits_week[visit_start:visit_end])

          # Walk reversed visits until day <= upper cutoff (mirrors Stan while loop)
          t_idx <- 1L
          while (t_idx <= n_visits && rev_days[t_idx] > max(0L, upper_sd[i_idx])) {
            t_idx <- t_idx + 1L
          }
          if (t_idx <= n_visits) {
            last_testing_visit_week[n, m, i] <- rev_weeks[t_idx]
            testing_end_idx[n, m, i]         <- visit_end - t_idx + 1L
          }
        }
      }
    }
  }

  list(
    first_testing_visit_week = first_testing_visit_week,
    last_testing_visit_week  = last_testing_visit_week,
    testing_start_idx        = testing_start_idx,
    testing_end_idx          = testing_end_idx
  )
}
