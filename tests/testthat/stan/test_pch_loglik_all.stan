functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}
data {
  int<lower=1> N;
  int<lower=1> MAX_T;
  int<lower=1> N_EXIT;
  array[N] int last_surv_week;
  array[N] int exit_event;
  array[N] int right_censored;
  array[N] int interval_censored;
  int ignore_interval_censoring;
  array[N_EXIT] matrix[N, MAX_T] log_cond_prob_surv;
  array[N] int start_from;
  array[N] int end_at;
}
generated quantities {
  // Full 8-arg version
  vector[N] lp_full = calc_pch_loglik(
    last_surv_week, exit_event, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv, start_from, end_at
  );
  // 7-arg: default end_at = max_all_t
  vector[N] lp_no_end = calc_pch_loglik(
    last_surv_week, exit_event, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv, start_from
  );
  // Single-type matrix overload (no exit_event, no start/end)
  vector[N] lp_single = calc_pch_loglik(
    last_surv_week, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv[1]
  );
  // pch_lpmf scalar sum — use sum(calc_pch_loglik(...)) to avoid pipe-syntax
  // restrictions on _lpmf in generated quantities assignments
  real lpmf_val = sum(calc_pch_loglik(
    last_surv_week, right_censored, interval_censored,
    ignore_interval_censoring, log_cond_prob_surv[1]
  ));
}
