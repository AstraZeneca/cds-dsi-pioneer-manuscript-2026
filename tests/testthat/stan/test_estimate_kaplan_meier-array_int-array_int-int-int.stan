functions {
  #include "../../../stan/pos.stan"
  #include "../../../stan/util.stan"
  #include "../../../stan/pfs_functions.stan"
}

data {
  int<lower=0> n_patients;
  array[n_patients] int<lower=0> event_time;
  array[n_patients] int<lower=0,upper=1> right_censored;
  int<lower=0> max_t;
  int<lower=0> pfs_offset;
}

generated quantities {
  vector[max_t + 1] km_survival;
  array[max_t + 1] int at_risk;
  array[max_t + 1] int n_right_censored;
  array[max_t + 1] int n_exited;
  
  (km_survival, at_risk, n_right_censored, n_exited) = estimate_kaplan_meier(event_time, right_censored, max_t, pfs_offset);
}
