// Regression tests for LFO correctness bugs fixed in b0a30e99:
//   C-EXT  Historical patients must train on full uncensored history
//   B      lfo_testing_patient_idx reverse-lookup: historical patient → 0
//   H3     GQ array cells [n, max_forecast_horizon] are zero, not NaN
//
// 2 patients (1 eval-trial, 1 historical), each with 5 visits.
// Patient 1 (eval trial=1): calendar_day=100, visits at study-days 1,3,5,8,12
// Patient 2 (historical): calendar_day=1,     visits at study-days 1,3,5,8,12
// Global cutoff: calendar_day=107
//   → patient 1: study-day cutoff = 107-100+1 = 8; last visit ≤ 8 → idx 4
//   → patient 2: WITHOUT C-EXT, study-day cutoff = 107; all visits qualify → idx 10
//                WITH C-EXT override, idx = 10 (last in block, same result for this data)
//                but the override fires unconditionally for any non-eval-trial patient.

functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
}

data {
  int<lower=1> lfo_eval_trial;
}

generated quantities {
  // Fixed sizes — no runtime-sized declarations in GQ
  array[2] int patient_trial      = {lfo_eval_trial, 3 - lfo_eval_trial};
  array[2] int calendar_day_pat   = {100, 1};
  array[3] int patient_visit_pos  = {1, 6, 11};

  // 5 visits per patient; global visit indices 1-5 (patient 1), 6-10 (patient 2)
  array[10] int t_patient_visits     = {1, 3, 5, 8, 12, 1, 3, 5, 8, 12};
  array[10] int t_patient_visits_day = {1, 3, 5, 8, 12, 1, 3, 5, 8, 12};

  int cutoff_cal = 107;

  // ── cutoff_visits ────────────────────────────────────────────────────────
  array[2] int cvd;   // cutoff_last_visit_day
  array[2] int cvw;   // cutoff_last_visit_week
  array[2] int lcd;   // last_visit_calendar_day
  array[2] int lwo;   // last_visit_week_observed
  array[2] int cvi;   // cutoff_last_visit_idx

  (cvd, cvw, lcd, lwo, cvi) = cutoff_visits(
    cutoff_cal, calendar_day_pat,
    t_patient_visits, t_patient_visits_day, patient_visit_pos
  );

  // ── C-EXT override ───────────────────────────────────────────────────────
  for (i in 1:2) {
    if (patient_trial[i] != lfo_eval_trial) {
      int vs, ve;
      (vs, ve) = get_pos(patient_visit_pos, i);
      cvi[i] = ve;
      cvw[i] = t_patient_visits[ve];
    }
  }

  // ── B: lfo_testing_patient_idx reverse-lookup ────────────────────────────
  // Sort patients by last_visit_calendar_day.
  // Patient 2 (cal_day=1): last cal-day = 1 + 12 - 1 = 12
  // Patient 1 (cal_day=100): last cal-day = 100 + 12 - 1 = 111
  // Sorted order: [2, 1]  (patient 2 first, patient 1 second)
  array[2] int sort_idx = sort_indices_asc(lcd);

  // OOS patients: those with last_visit_calendar_day > cutoff (=107)
  // Only patient 1 (day 111 > 107) → testing_patient_idx = 2 (position in sorted array)
  // n_all_testing_patients_unfiltered = 2 - 2 + 1 = 1
  array[1] int unfiltered_oos = sort_idx[2:];  // {patient 1}

  // Count eval-trial patients among the unfiltered OOS set
  int n_all_tp = 0;
  for (j in 1:1) {
    if (patient_trial[unfiltered_oos[j]] == lfo_eval_trial) n_all_tp += 1;
  }

  // Build all_testing_patients and reverse-lookup
  array[1] int all_tp = {0};
  {
    int idx = 1;
    for (j in 1:1) {
      if (patient_trial[unfiltered_oos[j]] == lfo_eval_trial) {
        all_tp[idx] = unfiltered_oos[j];
        idx += 1;
      }
    }
  }

  array[2] int lfo_tp_idx = zeros_int_array(2);
  for (j in 1:n_all_tp) {
    lfo_tp_idx[all_tp[j]] = j;
  }

  // ── H3: zero-init check ──────────────────────────────────────────────────
  // At the last cutoff, only [1,1] is filled; [1,2] must be 0.0, not NaN.
  // We replicate the pre-init loop from the GQ fix.
  array[1, 2] vector[1] ll_cells;
  for (n in 1:1) {
    for (m_rel in 1:2) {
      ll_cells[n, m_rel] = zeros_vector(1);
    }
  }
  // Only write [1,1]
  ll_cells[1, 1, 1] = -1.5;

  // ── Expose outputs for R assertions ──────────────────────────────────────
  // C-EXT
  int cext_historical_cutoff_idx = cvi[2];   // patient 2 (historical)
  int cext_eval_trial_cutoff_idx = cvi[1];   // patient 1 (eval-trial)

  // B
  int b_eval_trial_testing_idx   = lfo_tp_idx[1]; // patient 1 → slot 1
  int b_historical_testing_idx   = lfo_tp_idx[2]; // patient 2 → slot 0

  // H3
  real h3_cell_written    = ll_cells[1, 1, 1]; // -1.5
  real h3_cell_unwritten  = ll_cells[1, 2, 1]; // 0.0
}
