// tests/testthat/stan/test_multistate_ic.stan
// Tests for IC marginalization in multistate_lpmf
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  // Shared constants
  int n_patients;
  int max_t;
  int max_sojourn_t;

  // Patient state data
  array[n_patients] int final_state;
  array[n_patients] int time_01;
  array[n_patients] int time_02;
  array[n_patients] int time_12;
  array[n_patients] int time_03;
  array[n_patients] int time_32;
  array[n_patients] int censored_01;
  array[n_patients] int censored_02;
  array[n_patients] int censored_12;
  array[n_patients] int prog_deterministic;
  array[n_patients] int ms_ic_gap_01;

  // Feature flags
  int enable_03;

  // Visit arrays for visit-conditioned 0→3 dropout hazard
  int n_total_visits;
  array[n_total_visits] int t_patient_visits;
  array[n_patients + 1] int patient_visit_pos;

  // Constant log-conditional-survival (for tractable hand-computed expected values)
  real log_surv_val;   // single value — all hazards identical, all patients identical
}
transformed data {
  // Fill all hazard matrices with constant log_surv_val
  matrix[n_patients, max_t] lcs_01;
  matrix[n_patients, max_t] lcs_02;
  matrix[n_patients, max_t] lcs_03;
  matrix[n_patients, max_sojourn_t] lcs_12_s;
  matrix[n_patients, max_t] lcs_12_t;
  matrix[n_patients, max_t] lcs_32;
  for (i in 1:n_patients) {
    for (t in 1:max_t)         { lcs_01[i,t] = log_surv_val; lcs_02[i,t] = log_surv_val;
                                  lcs_03[i,t] = log_surv_val; lcs_12_t[i,t] = log_surv_val;
                                  lcs_32[i,t] = log_surv_val; }
    for (t in 1:max_sojourn_t) { lcs_12_s[i,t] = log_surv_val; }
  }
}
parameters { real dummy; }
model { dummy ~ normal(0, 1); }
generated quantities {
  real ll_ic = multistate_lpmf(
    final_state | 1, 1, 1, 1,   // enable_01, 02, 12, time_scale=semi-Markov
    enable_03, 0,                 // enable_03 (data), enable_32=0
    time_01, time_02, time_12,
    time_03, time_32,
    censored_01, censored_02, censored_12,
    rep_array(1, n_patients),    // censored_32 unused
    prog_deterministic,
    ms_ic_gap_01,
    t_patient_visits,
    patient_visit_pos,
    lcs_01, lcs_02, lcs_12_s, lcs_12_t, lcs_03, lcs_32,
    0, rep_vector(0.0, 0), 0.0    // enable_ms_visit_gated_01=0, no covariate
  );

  // Reference: same call with all gaps zeroed (current no-IC behavior)
  real ll_no_ic = multistate_lpmf(
    final_state | 1, 1, 1, 1,
    enable_03, 0,
    time_01, time_02, time_12,
    time_03, time_32,
    censored_01, censored_02, censored_12,
    rep_array(1, n_patients),
    prog_deterministic,
    rep_array(0, n_patients),    // gaps all zero -> no-IC path
    t_patient_visits,
    patient_visit_pos,
    lcs_01, lcs_02, lcs_12_s, lcs_12_t, lcs_03, lcs_32,
    0, rep_vector(0.0, 0), 0.0    // enable_ms_visit_gated_01=0, no covariate
  );
}
