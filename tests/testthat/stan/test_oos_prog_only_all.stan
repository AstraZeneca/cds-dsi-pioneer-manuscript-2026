// Tests for progression-only PFS in calculate_all_patients_endpoints_rng (Issue #92).
//
// Guards the LFO OOS RECIST confusion-matrix stamp: a direct-death (0->2) patient
// has sample_prog_right_censored=1 (progression-only PFS is censored, so the OOS
// stamp cannot write PD), while operational sample_right_censored=0 (death IS a
// PFS event). A 0->1 progression patient has sample_prog_right_censored=0.
//
// Run via test_stan_function() with fixed_param=TRUE, chains=1, iter_sampling=1.
// The R test asserts n_failures == 0.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "pfs.stanfunctions"
  #include "multistate.stanfunctions"
}

data {
  int<lower=0> dummy;  // Required non-empty data block; pass dummy = 0L from R
}

generated quantities {
  int n_failures = 0;

  // Build a 2-patient RNG-free scenario:
  //   P1: Direct death (0->2) without progression — operational PFS counts death
  //       as an event (sample_right_censored=0), but progression-only PFS must be
  //       censored (sample_prog_right_censored=1) so the OOS stamp cannot write PD.
  //   P2: Observed 0->1 progression — both operational and progression-only PFS
  //       report an uncensored event (sample_right_censored=0, sample_prog_right_censored=0).
  //
  // RNG elimination strategy (both patients):
  //   - Pass 0x0 ms_log_cond_surv → 0->1 sampling disabled for all patients.
  //   - P1: set right_censored[1]=1 (no observed 0->1), ms_censored_02[1]=0 (observed death)
  //         → sample_t01=max_all_t censored (line 1919-1920), sample_time_02=ms_time_02[1]
  //           deterministic (line 2020), cause=2 → derive_sample_pfs returns (exit_time, 0, ...)
  //         → operational sample_right_censored[1]=0 (death IS an event in operational PFS)
  //         BUT sample_prog_* captures sample_t01/sample_c01 BEFORE the death graft
  //         → sample_prog_pfs[1]=max_all_t, sample_prog_right_censored[1]=1 (censored).
  //   - P2: set right_censored[2]=0, pfs[2]=prog_week → sample_t01=pfs[2] deterministic
  //         (line 2049), sample_c01=0 → sample_prog_right_censored[2]=0.
  //   - forecast_size=0 for both patients → else branch (line 1861-1872) sets
  //     sample_target_pfs[j]=target_pfs[p] deterministically, no RNG.
  //   - Disable all unnecessary transitions: enable_ms_03=0, enable_ms_32=0, enable_ms_12=0.

  {
    int n = 2;
    int max_all_t = 50;

    // Shared arrays (0-sized or minimal where required)
    array[n] int patient_idx = { 1, 2 };
    array[4] int biomarker_category = { 1, 1, 1, 1 };  // 4 visits, SD or similar
    array[4] int rep_biomarker_category = { 1, 1, 1, 1 };  // 4 visits
    array[0] int forecast_obs_biomarker_category;  // Empty array
    array[n + 1] int forecast_obs_visits_pos = { 1, 1, 1 };  // No forecast obs visits (all ranges empty)
    int forecast_observation_interval = 6;
    matrix[0, max_all_t] ms_log_cond_surv;  // 0->1 disabled (0 patients worth of data)

    // Multistate enables: 02=ON (needed for P1 death), others OFF
    int enable_ms_02 = 1;
    int enable_ms_03 = 0;
    int enable_ms_32 = 0;
    int enable_ms_12 = 0;
    int ms_time_scale_12 = 0;
    int enable_ms_visit_gated_01 = 0;
    real tv_coef_01_val = 0.0;
    vector[0] forecast_obs_log_psa = zeros_vector(0);
    real median_log_psa_obs = 0.0;
    real iqr_log_psa_obs = 0.0;

    // 0->2 hazard (needed for P1 death); others 0 columns (disabled)
    matrix[n, max_all_t] log_cond_surv_02;
    for (i in 1:n) {
      for (t in 1:max_all_t) {
        log_cond_surv_02[i, t] = -0.05;  // Mild hazard
      }
    }
    matrix[n, 0] log_cond_surv_03;  // 0 columns = disabled
    matrix[n, 0] log_cond_surv_32;  // 0 columns = disabled
    matrix[n, 0] log_cond_surv_12_s;  // 0 columns = disabled
    matrix[n, 0] log_cond_surv_12_t;  // 0 columns = disabled

    // Patient 1: Direct death (0->2) at week 10, no progression
    // Patient 2: Observed 0->1 progression at week 12
    array[n] int pfs = { max_all_t, 12 };  // P1 censored for 0->1, P2 observed progression
    array[n] int right_censored = { 1, 0 };  // P1 censored for 0->1, P2 uncensored
    array[n] int target_pfs = { max_all_t, 12 };  // target-lesion PFS (mirrors pfs)
    array[n] int target_right_censored = { 1, 0 };
    array[n] int ms_right_censored = { 1, 0 };  // P1 censored for 0->1, P2 uncensored
    array[n] int ms_final_state = { 2, 1 };  // P1 died (state 2), P2 progressed (state 1)
    array[n] int ms_time_02 = { 10, 0 };  // P1 death at week 10
    array[n] int ms_censored_02 = { 0, 1 };  // P1 observed death, P2 no death
    array[n] int ms_censored_12 = { 1, 1 };  // Both censored for 1->2
    array[n] int ms_time_12 = { 0, 0 };
    array[n] int ms_os_event_12 = { 0, 0 };
    array[n] int ms_time_03 = { max_all_t, max_all_t };  // No dropout
    array[n] int ms_time_32 = { 0, 0 };
    array[n] int ms_censored_32 = { 1, 1 };

    // Visit position arrays: each patient has 2 visits (1 screening + 1 treatment)
    // to satisfy get_visit_pos validation: screening_visit_end >= visit_start
    // patient_visit_pos[i+1] - patient_visit_pos[i] = n_visits for patient i
    // P1: visits 1-2, P2: visits 3-4
    array[n + 1] int patient_visit_pos = { 1, 3, 5 };
    array[n + 1] int forecast_visits_pos = { 1, 1, 1 };  // No forecast visits
    array[n] int patient_last_obs_visit = { 6, 12 };  // Last observed visit weeks
    int last_predict_visit = max_all_t;
    array[4] int t_patient_visits = { 3, 6, 9, 12 };  // 2 visits per patient (screening + treatment)
    array[n] int n_patient_screening_visits = { 1, 1 };  // 1 screening visit each

    // Declare return variables
    array[n] int sample_target_pfs;
    array[n] int sample_target_right_censored;
    array[n] int spop_target_pfs;
    array[n] int spop_target_right_censored;
    array[n] int spop_target_obs_cens_pfs;
    array[n] int spop_target_obs_cens_right_censored;
    array[n] int sample_ms_pfs;
    array[n] int sample_ms_right_censored;
    array[n] int spop_ms_pfs;
    array[n] int spop_ms_right_censored;
    array[n] int sample_pfs;
    array[n] int sample_right_censored;
    array[n] int spop_pfs;
    array[n] int spop_right_censored;
    array[n] int sample_target_confirmed_response;
    array[n] int sample_target_unconfirmed_response;
    array[n] int spop_target_confirmed_response;
    array[n] int spop_target_unconfirmed_response;
    array[0] int forecast_target_pfs;
    array[0] int forecast_target_right_censored;
    array[n] int sample_os;
    array[n] int sample_os_censored;
    array[n] int spop_os;
    array[n] int spop_os_censored;
    array[n] int spop_is_dropout;
    array[n] int sample_is_dropout;
    array[n] int spop_dropout_week;
    array[n] int sample_prog_pfs;
    array[n] int sample_prog_right_censored;

    // Call calculate_all_patients_endpoints_rng
    (sample_target_pfs,
     sample_target_right_censored,
     spop_target_pfs,
     spop_target_right_censored,
     spop_target_obs_cens_pfs,
     spop_target_obs_cens_right_censored,
     sample_ms_pfs,
     sample_ms_right_censored,
     spop_ms_pfs,
     spop_ms_right_censored,
     sample_pfs,
     sample_right_censored,
     spop_pfs,
     spop_right_censored,
     sample_target_confirmed_response,
     sample_target_unconfirmed_response,
     spop_target_confirmed_response,
     spop_target_unconfirmed_response,
     forecast_target_pfs,
     forecast_target_right_censored,
     sample_os,
     sample_os_censored,
     spop_os,
     spop_os_censored,
     spop_is_dropout,
     sample_is_dropout,
     spop_dropout_week,
     sample_prog_pfs,
     sample_prog_right_censored) =
      calculate_all_patients_endpoints_rng(
        patient_idx,
        biomarker_category,
        rep_biomarker_category,
        forecast_obs_biomarker_category,
        forecast_obs_visits_pos,
        forecast_observation_interval,
        ms_log_cond_surv,
        enable_ms_02,
        enable_ms_03,
        enable_ms_32,
        enable_ms_12,
        ms_time_scale_12,
        log_cond_surv_02,
        log_cond_surv_03,
        log_cond_surv_32,
        log_cond_surv_12_s,
        log_cond_surv_12_t,
        pfs,
        right_censored,
        target_pfs,
        target_right_censored,
        ms_right_censored,
        ms_final_state,
        ms_time_02,
        ms_censored_02,
        ms_censored_12,
        ms_time_12,
        ms_os_event_12,
        ms_time_03,
        ms_time_32,
        ms_censored_32,
        patient_visit_pos,
        forecast_visits_pos,
        patient_last_obs_visit,
        last_predict_visit,
        t_patient_visits,
        max_all_t,
        n_patient_screening_visits,
        enable_ms_visit_gated_01,
        tv_coef_01_val,
        forecast_obs_log_psa,
        median_log_psa_obs,
        iqr_log_psa_obs
      );

    // Patient 1 checks: direct death (0->2) without progression
    // Operational PFS: death IS an event → sample_right_censored[1] == 0
    if (sample_right_censored[1] != 0) {
      print("FAIL P1 operational PFS: death is an event → sample_right_censored[1] should be 0, got ", sample_right_censored[1]);
      n_failures += 1;
    }
    // Progression-only PFS: no 0->1 event → sample_prog_right_censored[1] == 1 (censored)
    if (sample_prog_right_censored[1] != 1) {
      print("FAIL P1 progression-only PFS: no 0->1 → sample_prog_right_censored[1] should be 1, got ", sample_prog_right_censored[1]);
      n_failures += 1;
    }
    // Progression-only PFS time should be max_all_t (censored at horizon)
    if (sample_prog_pfs[1] != max_all_t) {
      print("FAIL P1 progression-only PFS time: expected max_all_t=", max_all_t, ", got ", sample_prog_pfs[1]);
      n_failures += 1;
    }

    // Patient 2 checks: observed 0->1 progression at week 12
    // Progression-only PFS: 0->1 event → sample_prog_right_censored[2] == 0
    if (sample_prog_right_censored[2] != 0) {
      print("FAIL P2 progression-only PFS: 0->1 event → sample_prog_right_censored[2] should be 0, got ", sample_prog_right_censored[2]);
      n_failures += 1;
    }
    // Progression-only PFS time should equal the progression week
    if (sample_prog_pfs[2] != 12) {
      print("FAIL P2 progression-only PFS time: expected 12, got ", sample_prog_pfs[2]);
      n_failures += 1;
    }
  }
}
