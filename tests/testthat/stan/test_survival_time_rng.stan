functions {
  #include "multistate.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_T;
  int<lower=1> N_DRAWS;
  array[N_CASES] int<lower=1> T;
  array[N_CASES] row_vector[MAX_T] log_cond_prob_surv;
  array[N_CASES] int<lower=0> obs_surv_time;
  array[N_CASES] int<lower=0,upper=1> right_censored;
  array[N_CASES] int<lower=0> interval_censored;
}

generated quantities {
  array[N_CASES, N_DRAWS] int sampled_time;
  array[N_CASES, N_DRAWS] int sampled_censored;
  for (i in 1:N_CASES) {
    for (j in 1:N_DRAWS) {
      tuple(int, int) res = survival_time_rng(log_cond_prob_surv[i, 1:T[i]], obs_surv_time[i], right_censored[i], interval_censored[i]);
      sampled_time[i, j] = res.1;
      sampled_censored[i, j] = res.2;
    }
  }
}
