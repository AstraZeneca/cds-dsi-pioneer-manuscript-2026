# tests/testthat/helper-multistate.R
#
# Pure R oracles for multistate log-likelihood functions.
# Mirrors calc_ms_single_transition_loglik and multistate_lpmf from
# stan/multistate.stanfunctions.

#' Oracle for calc_ms_single_transition_loglik
#'
#' Wraps the detection-week -> last-surviving-week convention shift before
#' computing the piecewise-constant hazard log-likelihood.
#'
#' @param event_time integer vector [n_patients] (detection-week convention)
#' @param censored   integer vector [n_patients] (1=censored, 0=event)
#' @param log_cond_surv matrix [n_patients x max_t] log conditional survival
#' @return numeric vector [n_patients] of per-patient log-likelihoods
r_calc_ms_stl <- function(event_time, censored, log_cond_surv) {
  n     <- length(event_time)
  max_t <- ncol(log_cond_surv)
  lp    <- numeric(n)

  for (i in seq_len(n)) {
    last_surv <- if (censored[i] == 1L) event_time[i] else event_time[i] - 1L

    # Known survival: weeks 1..last_surv
    if (last_surv >= 1L) {
      lp[i] <- sum(log_cond_surv[i, 1:last_surv])
    }

    # Hazard at event week (only for events, not censored)
    if (censored[i] == 0L) {
      # Stan calc_pch_loglik: ic_mix_lp[1] += log1m_exp(log_cond_prob_surv[1, i, last_surv + 1])
      ev_col <- last_surv + 1L
      if (ev_col >= 1L && ev_col <= max_t) {
        lp[i] <- lp[i] + log1p(-exp(log_cond_surv[i, ev_col]))
      }
    }
  }
  lp
}

#' Oracle for multistate_lpmf — 01-only fast path
#'
#' When enable_01=1, enable_02=0, enable_12=0 the Stan function delegates
#' directly to sum(calc_ms_single_transition_loglik(...)).
#'
#' @param time_01    integer vector [n_patients]
#' @param censored_01 integer vector [n_patients]
#' @param log_cond_surv_01 matrix [n_patients x max_t]
#' @return scalar log-probability
r_multistate_lpmf_01only <- function(time_01, censored_01, log_cond_surv_01) {
  sum(r_calc_ms_stl(time_01, censored_01, log_cond_surv_01))
}

#' Oracle for multistate_lpmf — state 0 censored patient with 01+02 enabled
#'
#' A right-censored (state-0) patient accumulates BOTH 01 and 02 survival
#' terms through the censoring time.
#'
#' @param t_cens integer scalar, censoring time
#' @param lcs_01 numeric vector [max_t] log conditional survival for 0->1
#' @param lcs_02 numeric vector [max_t] log conditional survival for 0->2
#' @return scalar log-probability contribution
r_multistate_state0_cens <- function(t_cens, lcs_01, lcs_02) {
  lp <- 0.0
  if (t_cens >= 1L) {
    lp <- lp + sum(lcs_01[1:t_cens]) + sum(lcs_02[1:t_cens])
  }
  lp
}
