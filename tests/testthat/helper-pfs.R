#' R reference implementation of Stan truncate_at_max_time.
#'
#' Mirrors: pfs.stanfunctions::truncate_at_max_time(pfs, right_censored, max_time)
#'
#' Semantics: if pfs[i] > max_time[i], truncate to max_time[i] and censor.
#' Boundary (pfs[i] == max_time[i]) is NOT truncated.
#'
#' @param pfs           Integer vector of PFS times.
#' @param right_censored Integer vector (0/1) of censoring flags.
#' @param max_time      Integer vector of per-patient cutoff times.
#' @return Named list with elements `pfs` and `right_censored`.
r_truncate_at_max_time <- function(pfs, right_censored, max_time) {
  n     <- length(pfs)
  o_pfs <- integer(n)
  o_rc  <- integer(n)
  for (i in seq_len(n)) {
    if (pfs[i] > max_time[i]) {
      o_pfs[i] <- max_time[i]
      o_rc[i]  <- 1L
    } else {
      o_pfs[i] <- pfs[i]
      o_rc[i]  <- right_censored[i]
    }
  }
  list(pfs = o_pfs, right_censored = o_rc)
}

#' R port of Stan find_first(all, what, n_succ).
#'
#' Returns the 1-based start index of the first run of n_succ consecutive
#' elements of `all` that are each present in `what`. Returns 0 if not found.
r_find_first <- function(all, what, n_succ) {
  n <- length(all)
  if (n < n_succ) return(0L)
  if (length(what) == 0L) stop("r_find_first: 'what' cannot be empty")
  if (n_succ <= 0L) stop("r_find_first: 'n_succ' must be >= 1")
  for (i in seq_len(n - n_succ + 1L)) {
    matches <- 0L
    for (j in seq_len(n_succ)) {
      if (all[i + j - 1L] %in% what) {
        matches <- matches + 1L
      } else {
        break
      }
    }
    if (matches == n_succ) return(i)
  }
  return(0L)
}

#' R port of Stan map_idx_to_week.
#'
#' Maps a 1-based combined (obs + forecast) index to an actual week number.
#' idx <= n_obs -> curr_visits[idx]; else -> forecast_time[idx - n_obs].
#' idx == 0 -> max_all_t (censoring sentinel).
r_map_idx_to_week <- function(idx, curr_visits, forecast_time, max_all_t) {
  n_obs <- length(curr_visits)
  if (idx == 0L) return(max_all_t)
  if (idx <= n_obs) return(curr_visits[idx])
  forecast_idx <- idx - n_obs
  if (forecast_idx < 1L || forecast_idx > length(forecast_time)) return(max_all_t)
  forecast_time[forecast_idx]
}

#' R reference implementation of Stan find_first_week.
#'
#' Mirrors: pfs.stanfunctions::find_first_week(arr, values, min_run_length,
#'                                              curr_visits, forecast_time, max_all_t)
#'
#' @param arr            Integer vector (RECIST codes or similar).
#' @param values         Integer vector of target values to search for.
#' @param min_run_length Minimum consecutive matches required.
#' @param curr_visits    Observed treatment visit weeks (may be integer(0)).
#' @param forecast_time  Forecast visit weeks.
#' @param max_all_t      Censoring sentinel returned when not found.
#' @return Named list with elements `week` (int) and `right_censored` (0/1 int).
r_find_first_week <- function(arr, values, min_run_length,
                              curr_visits, forecast_time, max_all_t) {
  idx <- r_find_first(arr, values, min_run_length)
  if (idx == 0L) {
    list(week = max_all_t, right_censored = 1L)
  } else {
    list(
      week          = r_map_idx_to_week(idx, curr_visits, forecast_time, max_all_t),
      right_censored = 0L
    )
  }
}

#' R reference implementation of Stan find_first_forecast_week.
#'
#' Delegates to r_find_first_week with empty curr_visits.
#'
#' @param arr            Integer vector.
#' @param values         Integer vector of target values.
#' @param min_run_length Minimum consecutive matches.
#' @param forecast_time  Forecast visit weeks.
#' @param max_all_t      Censoring sentinel.
#' @return Named list with `week` and `right_censored`.
r_find_first_forecast_week <- function(arr, values, min_run_length,
                                       forecast_time, max_all_t) {
  r_find_first_week(arr, values, min_run_length, integer(0), forecast_time, max_all_t)
}
