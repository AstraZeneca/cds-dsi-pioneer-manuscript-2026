
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

for (i in 1:n_patients) {
int visit_start, visit_end;
(visit_start, visit_end) = get_pos(patient_visit_pos, i); 

int start_idx = testing_start_idx[1, i];
// Only allocate OOS visits for patients who actually have a post-cutoff start.
// If start_idx == 0, the patient is not included at the first cutoff (no OOS window yet).
n_patient_testing_visits[i] = start_idx > 0 ? (visit_end - start_idx + 1) : 0;
}

testing_visit_pos = create_pos(n_patient_testing_visits);