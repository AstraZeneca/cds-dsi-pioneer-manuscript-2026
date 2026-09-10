# tests/testthat/helper-pch.R
#
# R oracle for calc_pch_loglik — mirrors the 8-arg Stan implementation exactly.
# log_cond_prob_surv_list: list of matrices, each [n_patients x max_t] (row = patient, col = week)

r_calc_pch_loglik <- function(
  last_surv_week, exit_event, right_censored, interval_censored,
  ignore_interval_censoring, log_cond_prob_surv_list,
  start_from, end_at
) {
  n_exit_types <- length(log_cond_prob_surv_list)
  n_patients   <- length(last_surv_week)
  lp <- numeric(n_patients)

  for (i in seq_len(n_patients)) {
    interval_pos <- max(0L, start_from[i])
    interval_end <- min(end_at[i], last_surv_week[i])

    eff_rc <- as.integer(
      right_censored[i] != 0L ||
      (end_at[i] < last_surv_week[i] + (1L - ignore_interval_censoring) * interval_censored[i] + 1L)
    )
    curr_ic <- if (ignore_interval_censoring != 0L || eff_rc != 0L) 0L else interval_censored[i]

    # Known-survival part: sum log_cond_surv over [interval_pos .. interval_end] (1-based cols)
    # Stan: sum(log_cond_prob_surv[k, i, interval_pos:interval_end])
    if (interval_pos <= interval_end) {
      for (k in seq_len(n_exit_types)) {
        lp[i] <- lp[i] + sum(log_cond_prob_surv_list[[k]][i, interval_pos:interval_end])
      }
    }

    # Interval-censoring mixture
    valid <- (interval_pos <= interval_end) ||
             (interval_end + curr_ic + 1L >= interval_pos)

    if (valid) {
      # Adjust interval_end and curr_ic (mirrors Stan's old_interval_end / old_interval_censored)
      old_interval_end <- interval_end
      interval_end     <- max(old_interval_end, interval_pos - 1L)
      curr_ic          <- max(0L, curr_ic - (interval_end - old_interval_end))

      mix <- numeric(curr_ic + 1L)

      for (c in 0:curr_ic) {
        if (c > 0L) {
          # Additional survival over (interval_end+1)..(interval_end+c)
          for (k in seq_len(n_exit_types)) {
            cols <- seq(interval_end + 1L, interval_end + c)
            if (length(cols) > 0L &&
                all(cols >= 1L) &&
                all(cols <= ncol(log_cond_prob_surv_list[[k]]))) {
              mix[c + 1L] <- mix[c + 1L] +
                sum(log_cond_prob_surv_list[[k]][i, cols])
            }
          }
        }

        if (eff_rc == 0L) {
          # Hazard at interval_end + c + 1
          # Stan: log1m_exp(log_cond_prob_surv[exit_event[i], i, interval_end + c + 1])
          col_ev <- interval_end + c + 1L
          if (col_ev >= 1L &&
              col_ev <= ncol(log_cond_prob_surv_list[[exit_event[i]]])) {
            mix[c + 1L] <- mix[c + 1L] +
              log1p(-exp(log_cond_prob_surv_list[[exit_event[i]]][i, col_ev]))
          }
        }
      }

      lp[i] <- lp[i] + if (curr_ic > 0L) {
        mat_max <- max(mix)
        mat_max + log(sum(exp(mix - mat_max))) - log(curr_ic + 1L)
      } else {
        mix[1L]
      }
    }
  }
  lp
}

# Convenience: single exit type, no IC, start=1, end=max_t
# log_cond_prob_surv: matrix [n_patients x max_t]
r_calc_pch_loglik_simple <- function(last_surv_week, right_censored, log_cond_prob_surv) {
  n_patients <- nrow(log_cond_prob_surv)
  max_t      <- ncol(log_cond_prob_surv)
  r_calc_pch_loglik(
    last_surv_week,
    rep(1L, n_patients),
    right_censored,
    rep(0L, n_patients),
    0L,
    list(log_cond_prob_surv),
    rep(1L, n_patients),
    rep(max_t, n_patients)
  )
}
