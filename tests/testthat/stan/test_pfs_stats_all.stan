functions {
  #include "multistate.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}
data {
  // calculate_log_marginal_exit_prob
  int<lower=1> T_lmep;
  row_vector[T_lmep] log_cond_surv;
  // survival_quantiles
  int<lower=1> N_sq;
  array[N_sq] int surv_time;
  int last_surv_time;
  int<lower=1> P_sq;
  vector[P_sq] sq_p;
  // km_quantiles + km_median + calc_km_pfs_n
  int<lower=1> T_km;
  vector[T_km] km_survival;
  int<lower=1> P_km;
  vector[P_km] km_p;
  real pfs_n_query;
}
generated quantities {
  row_vector[T_lmep] lmep = calculate_log_marginal_exit_prob(log_cond_surv);
  vector[P_sq] sq_quantiles;
  array[P_sq] int sq_cannot_calc;
  (sq_quantiles, sq_cannot_calc) = survival_quantiles(surv_time, last_surv_time, sq_p);
  real surv_med;
  int surv_med_flag;
  (surv_med, surv_med_flag) = survival_median(surv_time, last_surv_time);
  vector[P_km] km_q;
  array[P_km] int km_cannot_calc;
  (km_q, km_cannot_calc) = km_quantiles(km_survival, km_p);
  real km_med;
  int km_med_flag;
  (km_med, km_med_flag) = km_median(km_survival);
  real pfs_n = calc_km_pfs_n(km_survival, pfs_n_query);
}
