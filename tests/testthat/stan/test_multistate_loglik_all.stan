functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "multistate.stanfunctions"

  // Thin wrapper so we can call multistate_lpmf without the | pipe syntax
  // (which is rejected in generated quantities).
  real multistate_loglik(
    array[] int final_state,
    vector weight,
    int enable_01, int enable_02, int enable_12, int ms_time_scale_12,
    int enable_03, int enable_32,
    array[] int time_01, array[] int time_02, array[] int time_12,
    array[] int time_03, array[] int time_32,
    array[] int censored_01, array[] int censored_02, array[] int censored_12,
    array[] int censored_32,
    array[] int prog_deterministic,
    array[] int ms_ic_gap_01,
    array[] int t_patient_visits, array[] int patient_visit_pos,
    matrix log_cond_surv_01, matrix log_cond_surv_02,
    matrix log_cond_surv_12_s, matrix log_cond_surv_12_t,
    matrix log_cond_surv_03, matrix log_cond_surv_32
  ) {
    return multistate_lpmf(
      final_state |
      weight,
      enable_01, enable_02, enable_12, ms_time_scale_12,
      enable_03, enable_32,
      time_01, time_02, time_12, time_03, time_32,
      censored_01, censored_02, censored_12, censored_32,
      prog_deterministic,
      ms_ic_gap_01,
      t_patient_visits, patient_visit_pos,
      log_cond_surv_01, log_cond_surv_02,
      log_cond_surv_12_s, log_cond_surv_12_t,
      log_cond_surv_03, log_cond_surv_32
    );
  }
}
data {
  int<lower=1> N;
  int<lower=1> MAX_T;
  vector[N] weight;
  array[N] int event_time_01;
  array[N] int censored_01;
  matrix[N, MAX_T] log_cond_surv_01;
  // For full multistate test
  array[N] int final_state;
  array[N] int time_02;
  array[N] int censored_02;
  array[N] int time_12;
  array[N] int censored_12;
  matrix[N, MAX_T] log_cond_surv_02;
  matrix[N, MAX_T] log_cond_surv_12_s;
  matrix[N, MAX_T] log_cond_surv_12_t;
  array[N] int prog_deterministic;
  array[N] int ms_ic_gap_01;         // IC gap (0 = no interval censoring)
  // Visit arrays for 0→3 visit-conditioning (dummy when enable_03=0)
  int N_visits;
  array[N_visits] int t_patient_visits;
  array[N + 1] int patient_visit_pos;
  // Unused transitions (pass as sentinel)
  array[N] int time_03;
  array[N] int censored_32;
  matrix[N, MAX_T] log_cond_surv_03;
  matrix[N, MAX_T] log_cond_surv_32;
}
generated quantities {
  // Single-transition log-likelihood
  vector[N] st_llik = calc_ms_single_transition_loglik(
    event_time_01, censored_01, log_cond_surv_01
  );
  // multistate_lpmf: 01-only fast path
  real ms_lpmf_01only = multistate_loglik(
    final_state,
    weight,
    1, 0, 0, 0,
    0, 0,
    event_time_01, time_02, time_12, time_03, time_03,
    censored_01, censored_02, censored_12, censored_32,
    prog_deterministic,
    ms_ic_gap_01,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02,
    log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32
  );
  // multistate_lpmf: 01+02 (full illness-death, no 12)
  real ms_lpmf_01_02 = multistate_loglik(
    final_state,
    weight,
    1, 1, 0, 0,
    0, 0,
    event_time_01, time_02, time_12, time_03, time_03,
    censored_01, censored_02, censored_12, censored_32,
    prog_deterministic,
    ms_ic_gap_01,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02,
    log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32
  );
}
