// ============================================================================
// MS-Standalone LFO: Re-censor multistate events at calendar cutoff
// ============================================================================

// --- Step 1: Find each patient's last visit before the FIRST cutoff ---
array[n_patients] int cutoff_last_visit_day;
array[n_patients] int cutoff_last_visit_week;
array[n_patients] int last_visit_calendar_day;
array[n_patients] int last_visit_week_observed;
array[n_patients] int cutoff_last_visit_idx;

(cutoff_last_visit_day, cutoff_last_visit_week, last_visit_calendar_day,
 last_visit_week_observed, cutoff_last_visit_idx) =
  cutoff_visits(cutoff_calendar_day[1], calendar_day,
                t_patient_visits, t_patient_visits_day, patient_visit_pos);

// --- Step 2: Sort patients by last visit calendar day (for OOS indexing) ---
array[n_patients] int<lower=1, upper=n_patients> last_visit_calendar_day_sort_idx =
  sort_indices_asc(last_visit_calendar_day);

array[n_cutoffs] int<lower=1, upper=n_patients> testing_patient_idx =
  get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx],
                       cutoff_calendar_day);

print("testing_patient_idx = ", testing_patient_idx);
assert_ascending(testing_patient_idx);

int<lower=0, upper=n_patients> n_all_testing_patients_unfiltered =
  n_patients - testing_patient_idx[1] + 1;

// Filter testing patients to only those in the eval trial
int n_all_testing_patients = 0;
{
  array[n_all_testing_patients_unfiltered] int unfiltered =
    last_visit_calendar_day_sort_idx[testing_patient_idx[1]:];
  for (j in 1:n_all_testing_patients_unfiltered) {
    if (patient_trial[unfiltered[j]] == lfo_eval_trial) n_all_testing_patients += 1;
  }
}
array[n_all_testing_patients] int<lower=1, upper=n_patients> all_testing_patients;
{
  array[n_all_testing_patients_unfiltered] int unfiltered =
    last_visit_calendar_day_sort_idx[testing_patient_idx[1]:];
  int idx = 1;
  for (j in 1:n_all_testing_patients_unfiltered) {
    if (patient_trial[unfiltered[j]] == lfo_eval_trial) {
      all_testing_patients[idx] = unfiltered[j];
      idx += 1;
    }
  }
}

// Reverse lookup: patient i → position in all_testing_patients (0 = not in eval trial)
array[n_patients] int<lower=0, upper=n_all_testing_patients> lfo_testing_patient_idx =
  zeros_int_array(n_patients);
for (j in 1:n_all_testing_patients) {
  lfo_testing_patient_idx[all_testing_patients[j]] = j;
}

int<lower=0, upper=n_patients> n_training_patients = testing_patient_idx[1] - 1;

print("n_all_testing_patients (trial ", lfo_eval_trial, ") = ", n_all_testing_patients);
print("n_training_patients = ", n_training_patients);

// --- Step 3: Compute testing visit bounds for GQ evaluation windows ---
array[n_cutoffs, n_patients] int<lower=0> oos_patient_first_testing_visit_week;
array[n_cutoffs, n_cutoffs, n_patients] int oos_patient_last_testing_visit_week;
array[n_cutoffs, n_patients] int testing_start_idx;
array[n_cutoffs, n_cutoffs, n_patients] int testing_end_idx;

(oos_patient_first_testing_visit_week, oos_patient_last_testing_visit_week,
 testing_start_idx, testing_end_idx) =
  get_testing_visit_week_bounds(testing_patient_idx, last_visit_calendar_day_sort_idx,
                                cutoff_calendar_day, calendar_day,
                                t_patient_visits, t_patient_visits_day, patient_visit_pos);

// --- Step 4: Mark which patients are enrolled at the cutoff ---
// The actual lfo_likelihood_weight is computed in transformed parameters
// (after propensity likelihood_weight is available).
array[n_patients] int lfo_patient_enrolled;
for (i in 1:n_patients) {
  lfo_patient_enrolled[i] = cutoff_last_visit_idx[i] > 0 ? 1 : 0;
}

// --- Step 5: Compute per-patient calendar cutoff in week scale ---
// Death, dropout, and off-trial death are registry-exact events that do not
// require a visit to be observed, so they must be censored at the calendar
// cutoff — not the last visit before it. Visit-gated events (progression,
// state-0 follow-up) continue to use cutoff_last_visit_week.
array[n_patients] int lfo_cutoff_cal_week;
for (i in 1:n_patients) {
  int days_since_enroll = cutoff_calendar_day[1] - calendar_day[i] + 1;
  lfo_cutoff_cal_week[i] = days_since_enroll > 0 ? (days_since_enroll - 1) %/% 7 + 1 : 0;
}

// --- Step 6: Re-censor multistate events at the first cutoff ---
// Calls the shared recensor_ms_at_cutoff() function (tested in test_recensor_ms_all.stan).
array[n_patients] int lfo_ms_final_state;
array[n_patients] int lfo_ms_time_01;
array[n_patients] int lfo_ms_censored_01;
array[n_patients] int lfo_ms_time_02;
array[n_patients] int lfo_ms_time_12;
array[n_patients] int lfo_ms_time_03;
array[n_patients] int lfo_ms_time_32;
array[n_patients] int lfo_interval_censored;
array[n_patients] int lfo_ms_prog_deterministic;
array[n_patients] int lfo_ms_ic_gap_01;

(lfo_ms_final_state, lfo_ms_time_01, lfo_ms_censored_01,
 lfo_ms_time_02, lfo_ms_time_12, lfo_ms_time_03, lfo_ms_time_32,
 lfo_interval_censored, lfo_ms_prog_deterministic, lfo_ms_ic_gap_01) =
  recensor_ms_at_cutoff(
    ms_final_state, ms_time_01, ms_censored_01,
    ms_time_02, ms_time_12, ms_time_03, ms_time_32,
    ms_os_event_12, interval_censored, ms_prog_deterministic,
    cutoff_last_visit_week, lfo_cutoff_cal_week);
