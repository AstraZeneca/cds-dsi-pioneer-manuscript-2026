functions {
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
}
data {
  int N;
  int MAX_T;
  // Raw (un-censored) multistate fields
  array[N] int ms_final_state;
  array[N] int ms_time_01;
  array[N] int ms_censored_01;
  array[N] int ms_time_02;
  array[N] int ms_time_12;
  array[N] int ms_time_03;
  array[N] int ms_time_32;
  array[N] int ms_os_event_12;
  array[N] int interval_censored;
  array[N] int ms_prog_deterministic;
  array[N] int cutoff_visit_week;   // per-patient last visit at/before cutoff
  array[N] int cutoff_cal_week;     // per-patient calendar cutoff (week)
  // Visit arrays (enable_03 path uses these; keep simple, all weeks distinct)
  array[N + 1] int patient_visit_pos;   // length N+1
  array[size(patient_visit_pos) > 0 ? patient_visit_pos[N + 1] - 1 : 0] int t_patient_visits;
  // Conditional survival matrices (already -exp transformed log conditional survival)
  matrix[N, MAX_T] log_cond_surv_01;
  matrix[N, MAX_T] log_cond_surv_02;
  matrix[N, MAX_T] log_cond_surv_12_s;
  matrix[N, MAX_T] log_cond_surv_12_t;
  matrix[N, MAX_T] log_cond_surv_03;
  matrix[N, MAX_T] log_cond_surv_32;
  // Flags
  int enable_01; int enable_02; int enable_12; int ms_time_scale_12;
  int enable_03; int enable_32;
  int enable_visit_gated_01;   // 0 = continuous path, 1 = production visit-gated path
}
parameters { real dummy; }
model { dummy ~ std_normal(); }
generated quantities {
  // Re-censor at cutoff
  array[N] int f; array[N] int t01; array[N] int c01;
  array[N] int t02; array[N] int t12; array[N] int t03; array[N] int t32;
  array[N] int ic; array[N] int pd; array[N] int gap01;
  (f, t01, c01, t02, t12, t03, t32, ic, pd, gap01) =
    recensor_ms_at_cutoff(
      ms_final_state, ms_time_01, ms_censored_01,
      ms_time_02, ms_time_12, ms_time_03, ms_time_32,
      ms_os_event_12, interval_censored, ms_prog_deterministic,
      cutoff_visit_week, cutoff_cal_week);

  // LFO likelihood over re-censored data
  real lfo_ll = multistate_lpmf(
    f | ones_vector(N),
    enable_01, enable_02, enable_12, ms_time_scale_12, enable_03, enable_32,
    t01, t02, t12, t03, t32, c01, pd, gap01,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02, log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32,
    // enable_ms_visit_gated_01: threaded from data so the same scenario can be
    // run continuous-OFF (01/02/12 subset) or gated-ON (= production config,
    // enable_ms_visit_gated_01=1). Both lfo_ll and full_ll must use the SAME
    // value for the parity assertion to be meaningful.
    enable_visit_gated_01);

  // Full-model likelihood over the ORIGINAL (un-censored) data.
  // Parity edge verified: the prog_deterministic=1 patient (P3) yields
  // ic_gap_01 = 0 from BOTH recensor_ms_at_cutoff (gap01 above) and this
  // independent gap01_full derivation, so the two gap rules agree — this is
  // the spot where recensor and the full model are most likely to silently
  // diverge.
  array[N] int gap01_full;
  for (i in 1:N) gap01_full[i] = (ms_censored_01[i] || ms_prog_deterministic[i]) ? 0 : interval_censored[i] + 1;
  real full_ll = multistate_lpmf(
    ms_final_state | ones_vector(N),
    enable_01, enable_02, enable_12, ms_time_scale_12, enable_03, enable_32,
    ms_time_01, ms_time_02, ms_time_12, ms_time_03, ms_time_32,
    ms_censored_01, ms_prog_deterministic, gap01_full,
    t_patient_visits, patient_visit_pos,
    log_cond_surv_01, log_cond_surv_02, log_cond_surv_12_s, log_cond_surv_12_t,
    log_cond_surv_03, log_cond_surv_32,
    enable_visit_gated_01);  // same flag as lfo_ll above
}
