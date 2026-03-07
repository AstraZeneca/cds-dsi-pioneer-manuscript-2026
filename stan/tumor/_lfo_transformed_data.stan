
// --- LFO CV specific (visit-based) ---
array[n_patients] int cutoff_last_visit_day;
array[n_patients] int cutoff_last_visit_week;
array[n_patients] int last_visit_calendar_day;
array[n_patients] int last_visit_week_observed;

// Precompute for each patient the last visit index before or at cutoff_last_visit_week
array[n_patients] int cutoff_last_visit_idx;

(cutoff_last_visit_day, cutoff_last_visit_week, last_visit_calendar_day, last_visit_week_observed, cutoff_last_visit_idx) = cutoff_visits(
cutoff_calendar_day[1], calendar_day, t_patient_visits, t_patient_visits_day, patient_visit_pos
);

// This is an array of patient IDs (sorted by last visit calendar day)
array[n_patients] int<lower = 1, upper = n_patients> last_visit_calendar_day_sort_idx = sort_indices_asc(last_visit_calendar_day);

// Per cutoff, which index in the above sorted list of patient IDs, identifying the first patient to be included in the out-of-sample testing set
array[n_cutoffs] int<lower = 1, upper = n_patients> testing_patient_idx = 
get_oos_patients_idx(last_visit_calendar_day[last_visit_calendar_day_sort_idx], cutoff_calendar_day);

print("testing_patient_idx = ", testing_patient_idx);

assert_ascending(testing_patient_idx);

int<lower = 0, upper = n_patients> n_all_testing_patients = n_patients - testing_patient_idx[1] + 1;

// These are the patient IDs of all patients that are included in the out-of-sample testing set, sorted by last visit calendar day
array[n_all_testing_patients] int<lower = 1, upper = n_patients> all_testing_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[1]:];

print("n_all_testing_patients = ", n_all_testing_patients);

// Training patients are those observed before the first cutoff
int<lower = 0, upper = n_patients> n_training_patients = testing_patient_idx[1] - 1;

// These are the patient IDs of all patients in the training set, sorted by last visit calendar day
array[n_training_patients] int<lower = 1, upper = n_patients> training_patients = 
  n_training_patients > 0 ? last_visit_calendar_day_sort_idx[1:n_training_patients] : rep_array(0, 0);

print("n_training_patients = ", n_training_patients);

array[n_cutoffs, n_patients] int<lower = 0> oos_patient_first_testing_visit_week;
array[n_cutoffs, n_cutoffs, n_patients] int oos_patient_last_testing_visit_week;
array[n_cutoffs, n_patients] int testing_start_idx;
// testing_end_idx[n, m, i] stores the last inclusive visit index for patient i when evaluating
// a window that starts at cutoff n and ends just BEFORE cutoff m (i.e. m is the next cutoff).
// Therefore, for an evaluation horizon ending at cutoff m (with m >= n), we look up
// testing_end_idx[n, m + 1, i] unless m == n_cutoffs, in which case we fall back to the patient's
// final visit. This "shift by +1" in the second dimension lets us treat the final horizon uniformly
// without allocating an out-of-range m+1 cell.
array[n_cutoffs, n_cutoffs, n_patients] int testing_end_idx;

(oos_patient_first_testing_visit_week, oos_patient_last_testing_visit_week, testing_start_idx, testing_end_idx) = get_testing_visit_week_bounds(
testing_patient_idx, last_visit_calendar_day_sort_idx, cutoff_calendar_day, calendar_day, t_patient_visits, t_patient_visits_day, patient_visit_pos
);

for (n in 1:n_cutoffs) {
int n_curr_patients = n_patients - testing_patient_idx[n] + 1; // How many patients after the current patient index
array[n_curr_patients] int curr_patients = last_visit_calendar_day_sort_idx[testing_patient_idx[n]:]; // Who are these patients
array[n_cutoffs] int m_size = rep_array(-1, n_cutoffs);

for (m in 1:n_cutoffs) {
    if (m >= n) {
    m_size[m] = 0;

    for (i_idx in 1:n_curr_patients) {
        // Note: i is the original patient ID (1-based index from input data), not a sort position.
        // curr_patients contains original patient IDs that were reordered by sorting on last_visit_calendar_day
        int i = curr_patients[i_idx];

        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i); 

        int start_idx = testing_start_idx[n, i];
        int end_idx = m < n_cutoffs ? testing_end_idx[n, m + 1, i] : visit_end;

        if (start_idx > 0 && end_idx >= start_idx) {
        m_size[m] += 1;    
        }
    }
    }
}

print("Number of future visits:");
print(n,": ", m_size);
}

array[n_patients + 1] int<lower = 1> testing_visit_pos = zeros_int_array(n_patients + 1); 
array[n_patients] int n_patient_testing_visits;

// First, identify which patients have any observations at the cutoff
// (cutoff_last_visit_idx[i] > 0 means patient i has at least one visit before/at cutoff)
array[n_patients] int cutoff_observed_mask;
int n_cutoff_observed_patients = 0;

for (i in 1:n_patients) {
  int visit_start, visit_end;
  (visit_start, visit_end) = get_pos(patient_visit_pos, i); 
  
  int start_idx = testing_start_idx[1, i];
  
  // Identify patients observed at cutoff
  if (cutoff_last_visit_idx[i] > 0) {
    cutoff_observed_mask[i] = 1;
    n_cutoff_observed_patients += 1;
  } else {
    cutoff_observed_mask[i] = 0;
  }
  
  // Only allocate OOS visits for patients who:
  // 1) Have post-cutoff visits (start_idx > 0), AND
  // 2) Were observed before/at the cutoff (cutoff_observed_mask[i] == 1)
  // This matches the condition in sf-ssls-lfo.stan where predictions are generated
  n_patient_testing_visits[i] = (start_idx > 0 && cutoff_observed_mask[i] == 1) ? (visit_end - start_idx + 1) : 0;
}

print("n_cutoff_observed_patients = ", n_cutoff_observed_patients);

testing_visit_pos = create_pos(n_patient_testing_visits);


// Create a compact array of patient IDs who were observed at cutoff
// and a mapping from original patient ID to compact index
array[n_cutoff_observed_patients] int cutoff_observed_patients;
array[n_patients] int patient_to_cutoff_idx;
(cutoff_observed_patients, patient_to_cutoff_idx) = create_compact_patient_mapping(
  cutoff_observed_mask, n_cutoff_observed_patients
);

// Create cutoff-censored data ONLY for observed patients (compact arrays)
array[n_cutoff_observed_patients] int cutoff_pfs;
array[n_cutoff_observed_patients] int cutoff_right_censored;
array[n_cutoff_observed_patients] int cutoff_interval_censored;

for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];  // Original patient ID
  
  // If patient's observed PFS is after the cutoff, censor them at cutoff
  // cutoff_last_visit_week[i] is the last week observed for patient i at cutoff
  if (pfs[i] > cutoff_last_visit_week[i]) {
    cutoff_pfs[obs_idx] = cutoff_last_visit_week[i];
    cutoff_right_censored[obs_idx] = 1;  // Censored at cutoff
    cutoff_interval_censored[obs_idx] = 0;  // No interval censoring for cutoff-censored patients
  } else {
    // Event occurred before cutoff, use actual observed data
    cutoff_pfs[obs_idx] = pfs[i];
    cutoff_right_censored[obs_idx] = right_censored[i];
    cutoff_interval_censored[obs_idx] = interval_censored[i];
  }
}

// Create cutoff-censored multistate 0→1 data
array[n_cutoff_observed_patients] int cutoff_ms_time_01;
array[n_cutoff_observed_patients] int cutoff_ms_censored_01;

for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];  // Original patient ID

  // If patient's multistate 0→1 time is after the cutoff, censor them at cutoff
  if (ms_time_01[i] > cutoff_last_visit_week[i]) {
    cutoff_ms_time_01[obs_idx] = cutoff_last_visit_week[i];
    cutoff_ms_censored_01[obs_idx] = 1;  // Censored at cutoff
  } else {
    // Event occurred before cutoff, use actual observed data
    cutoff_ms_time_01[obs_idx] = ms_time_01[i];
    cutoff_ms_censored_01[obs_idx] = ms_censored_01[i];
  }
}

// Create cutoff-censored target lesion PFS data
array[n_cutoff_observed_patients] int cutoff_target_pfs;
array[n_cutoff_observed_patients] int cutoff_target_right_censored;

for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];  // Original patient ID
  
  // If patient's target PFS is after the cutoff, censor them at cutoff
  if (target_pfs[i] > cutoff_last_visit_week[i]) {
    cutoff_target_pfs[obs_idx] = cutoff_last_visit_week[i];
    cutoff_target_right_censored[obs_idx] = 1;  // Censored at cutoff
  } else {
    // Event occurred before cutoff, use actual observed data
    cutoff_target_pfs[obs_idx] = target_pfs[i];
    cutoff_target_right_censored[obs_idx] = target_right_censored[i];
  }
}

// Create compact patient state data for cutoff-observed patients
// First, count total visits for observed patients at cutoff
int n_cutoff_visits = 0;
array[n_cutoff_observed_patients] int cutoff_n_patient_visits;

for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];
  int last_cutoff_idx = cutoff_last_visit_idx[i];
  int start_idx = patient_visit_pos[i];
  
  // Count visits up to and including the cutoff
  int n_visits_at_cutoff = last_cutoff_idx > 0 ? (last_cutoff_idx - start_idx + 1) : 0;
  cutoff_n_patient_visits[obs_idx] = n_visits_at_cutoff;
  n_cutoff_visits += n_visits_at_cutoff;
}

// Compute n_cutoff_visits_m1 (excluding first visit per patient)
int n_cutoff_visits_m1 = n_cutoff_visits > n_cutoff_observed_patients ? n_cutoff_visits - n_cutoff_observed_patients : 0;

// Create compact position array for observed patients' visits
array[n_cutoff_observed_patients + 1] int cutoff_patient_visit_pos = create_pos(cutoff_n_patient_visits);

// Create compact RECIST array (only visits for observed patients before cutoff)
array[n_cutoff_visits] int cutoff_recist;

// Create compact arrays for all visit-level and patient-level data needed for state generation
// These will be used to call generate_all_patients_states_with_means_rng with only cutoff-observed data
array[n_cutoff_observed_patients] int cutoff_n_patient_screening_visits;
array[n_cutoff_observed_patients] int cutoff_patient_last_obs_visit;
array[n_cutoff_visits] int cutoff_t_patient_visits;

{
  int cutoff_visit_idx = 1;
  for (obs_idx in 1:n_cutoff_observed_patients) {
    int i = cutoff_observed_patients[obs_idx];
    int start_idx = patient_visit_pos[i];
    int last_cutoff_idx = cutoff_last_visit_idx[i];
    
    // Copy patient-level metadata
    cutoff_n_patient_screening_visits[obs_idx] = n_patient_screening_visits[i];
    cutoff_patient_last_obs_visit[obs_idx] = cutoff_n_patient_visits[obs_idx]; // Last visit at cutoff
    
    // Copy RECIST values and visit times for visits before/at cutoff
    for (v in start_idx:last_cutoff_idx) {
      cutoff_recist[cutoff_visit_idx] = recist[v];
      cutoff_t_patient_visits[cutoff_visit_idx] = t_patient_visits[v];
      cutoff_visit_idx += 1;
    }
  }
}

// Create cutoff-specific trial patient groupings (using compact patient indices)
array[n_trials + 1] int cutoff_trial_patient_pos = create_compact_group_pos(
  cutoff_observed_patients, patient_trial, n_trials
);

// Create cutoff-specific conditional group patient groupings (only include observed patients)
// Note: cond_group contains patient IDs (potentially with duplicates for multiple groups)
// First compute the size by counting observed entries
int n_cutoff_cond_group_entries = 0;
for (g_idx in 1:size(cond_group)) {
  if (cutoff_observed_mask[cond_group[g_idx]] == 1) {
    n_cutoff_cond_group_entries += 1;
  }
}

array[n_cutoff_cond_group_entries] int cutoff_cond_group;
array[n_cond_group + 1] int cutoff_cond_group_pos;
(cutoff_cond_group, cutoff_cond_group_pos) = remap_group_to_compact(
  cond_group, cond_group_pos, cutoff_observed_mask, patient_to_cutoff_idx, n_cond_group
);

// Forecast visits for patients right-censored at cutoff
// (either truly censored or censored because next visit is after cutoff)
array[n_cutoff_observed_patients] int cutoff_n_patient_forecast_visits;
int n_cutoff_right_censored_patients = 0;

for (obs_idx in 1:n_cutoff_observed_patients) {
  int i = cutoff_observed_patients[obs_idx];
  // cutoff_patient_last_obs_visit is already computed above in terms of visit count
  // Now compute the actual time of the last observation
  int last_obs_time = cutoff_patient_last_obs_visit[obs_idx] > 0 ? cutoff_t_patient_visits[cutoff_patient_visit_pos[obs_idx + 1] - 1] : 0;
  cutoff_n_patient_forecast_visits[obs_idx] = max_all_t - last_obs_time;
  
  // Count right-censored patients (either truly censored or their event/next visit is after cutoff)
  if (cutoff_right_censored[obs_idx] == 1) {
    n_cutoff_right_censored_patients += 1;
  }
}

int n_cutoff_total_forecast_visits = sum(cutoff_n_patient_forecast_visits);
array[n_cutoff_observed_patients + 1] int cutoff_forecast_visits_pos = create_pos(cutoff_n_patient_forecast_visits);

// Create mapping from cutoff visits to original state indices for subsetting in generated quantities
// Include ALL visits (including first visit per patient) - the function expects full visit states
array[n_cutoff_visits] int cutoff_state_indices;
{
  int cutoff_idx = 1;
  for (obs_idx in 1:n_cutoff_observed_patients) {
    int i = cutoff_observed_patients[obs_idx];
    int start_idx = patient_visit_pos[i];
    int last_cutoff_idx = cutoff_last_visit_idx[i];
    
    // Map each cutoff visit (INCLUDING first visit) to its corresponding state row index
    for (v in start_idx:last_cutoff_idx) {  // Include first visit
      cutoff_state_indices[cutoff_idx] = v;  // states uses absolute visit indexing
      cutoff_idx += 1;
    }
  }
}

// Create m1 position array (needed by the function for process noise indexing)
// This is derived directly from visit counts using create_pos with -1 increment
array[n_cutoff_observed_patients + 1] int cutoff_patient_visit_m1_pos = create_pos(cutoff_n_patient_visits, -1);

// Note: These patient-level parameters come from transformed parameters and will be subsetted in generated quantities
// where we have access to the parameter values using fancy indexing:
// - patient_log_decrease_rate[cutoff_observed_patients]  (patient-level)
// - patient_log_growth_rate[cutoff_observed_patients]    (patient-level)
// - sum_tumor_size[cutoff_state_indices]                 (visit-level!)
// - states[cutoff_state_indices, ] to get cutoff states (ALL visits, including first)

