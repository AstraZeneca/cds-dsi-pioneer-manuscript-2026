# tests/testthat/helper-pfs-stats.R
#
# Pure R oracles for pfs.stanfunctions statistical utilities.

r_calculate_log_marginal_exit_prob <- function(log_cond_prob_surv) {
  T   <- length(log_cond_prob_surv)
  out <- numeric(T)
  for (t in seq_len(T)) {
    out[t] <- log1p(-exp(log_cond_prob_surv[t]))
    if (t > 1) out[t] <- out[t] + sum(log_cond_prob_surv[seq_len(t - 1L)])
  }
  out
}

r_survival_quantiles <- function(surv_time, last_surv_time, p) {
  N         <- length(surv_time)
  sorted    <- sort(surv_time)
  quantiles    <- numeric(length(p))
  cannot_calc  <- integer(length(p))
  for (j in seq_along(p)) {
    k <- 1L
    while (p[j] >= k * (1.0 / N)) k <- k + 1L
    if (sorted[max(k - 1L, 1L)] > last_surv_time) {
      quantiles[j] <- 0.0; cannot_calc[j] <- 1L
    } else {
      pos <- p[j] * (N - 1) + 1
      d   <- pos - (k - 1)
      quantiles[j] <- sorted[max(k - 1L, 1L)] +
                      d * (sorted[k] - sorted[max(k - 1L, 1L)])
      cannot_calc[j] <- 0L
    }
  }
  list(quantiles = quantiles, cannot_calculate = cannot_calc)
}

r_km_quantiles <- function(km_survival, p) {
  T    <- length(km_survival)
  P    <- length(p)
  quantiles   <- rep(0.0, P)
  cannot_calc <- rep(1L, P)
  p_idx <- order(p, decreasing = TRUE)
  curr  <- 1L
  while (curr <= P && p[p_idx[curr]] > km_survival[1]) {
    quantiles[p_idx[curr]]   <- 0.0
    cannot_calc[p_idx[curr]] <- 0L
    curr <- curr + 1L
  }
  for (t in 2:T) {
    while (curr <= P &&
           km_survival[t]     <= p[p_idx[curr]] &&
           km_survival[t - 1] >  p[p_idx[curr]]) {
      idx    <- p_idx[curr]
      weight <- (km_survival[t - 1] - p[idx]) /
                (km_survival[t - 1] - km_survival[t])
      quantiles[idx]   <- (t - 2) + weight   # 0-based weeks
      cannot_calc[idx] <- 0L
      curr <- curr + 1L
    }
    if (curr > P) break
  }
  list(quantiles = quantiles, cannot_calculate = cannot_calc)
}

r_calc_km_pfs_n <- function(km_survival, n) {
  max_t <- length(km_survival) - 1L
  if (n <= 0)      return(km_survival[1L])
  if (n >= max_t)  return(km_survival[max_t + 1L])
  t_floor <- floor(n)
  weight  <- n - t_floor
  (1 - weight) * km_survival[t_floor + 1L] + weight * km_survival[t_floor + 2L]
}
