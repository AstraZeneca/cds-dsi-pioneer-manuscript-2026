// Multi-level hierarchy: computed totals and position arrays
int n_total_groups = sum(n_groups_per_level);
array[n_levels + 1] int level_pos = create_pos(n_groups_per_level);

// Validate multi-level hierarchy structure
for (p in 1:n_patients) {
  for (lv in 1:n_levels) {
    if (patient_level_groups[p, lv] < 1 ||
        patient_level_groups[p, lv] > n_groups_per_level[lv]) {
      fatal_error("Patient ", p, " has invalid group assignment at level ", lv,
                  ": got ", patient_level_groups[p, lv],
                  ", expected 1-", n_groups_per_level[lv]);
    }
  }
  // Verify patient level is identity mapping
  if (patient_level_groups[p, n_levels] != p) {
    fatal_error("Patient ", p, " must have identity mapping at patient level (level ",
                n_levels, "), got ", patient_level_groups[p, n_levels]);
  }
}

// Validate patient level has correct group count
if (n_groups_per_level[n_levels] != n_patients) {
  fatal_error("Patient level must have n_groups = n_patients, got ",
              n_groups_per_level[n_levels], " vs ", n_patients);
}

// Validate backward compatibility: patient_level_groups[,1] must match patient_trial
for (p in 1:n_patients) {
  if (patient_level_groups[p, 1] != patient_trial[p]) {
    fatal_error("patient_level_groups[", p, ", 1] = ", patient_level_groups[p, 1],
                " must match patient_trial[", p, "] = ", patient_trial[p]);
  }
}

print("Multi-level hierarchy validated:");
print("  n_levels = ", n_levels);
print("  n_groups_per_level = ", n_groups_per_level);

// Starting position of patients for each trial in a flattened patient array
// Diagram for trial_patient_pos:
// [1, 40, 70, 120, ...]
//  ^  ^   ^   ^
//  |  |   |   |
//  |  |   |   Start of patients in trial 4
//  |  |   Start of patients in trial 3
//  |  Start of patients in trial 2
//  Start of patients in trial 1
array[n_trials + 1] int<lower = 1, upper = n_patients + 1> trial_patient_pos;
{
  array[n_trials] int n_trial_patients = rep_array(0, n_trials);
  for (i in 1:n_patients) {
    n_trial_patients[patient_trial[i]] += 1;
  }
  print("n_trial_patients = ", n_trial_patients);
  trial_patient_pos = create_pos(n_trial_patients);
}

array[n_patients + 1] int<lower = 1, upper = sum(n_patient_visits) + 1> patient_visit_pos = create_pos(n_patient_visits);

real delta = 1e-5; // Small value used for GP modeling to avoid numerical issues

int max_all_t = max(max(t_patient_visits) + 1, extend_max_all_t); // Latest measurement time or extended time, whichever is greater
int<lower = 0> max_t_width = max_all_t - min(t_patient_visits) + 1;

print("max(t_patient_visits) = ", max(t_patient_visits));
print("max_all_t = ", max_all_t);
print("max_t_width = ", max_t_width);

array[sum(n_patient_visits)] int<lower = 1> t_patient_visit_idx; // Index of each patient visit relative to the first visit for each patient 

array[n_patients] int<lower = 0> n_patient_screening_visits = zeros_int_array(n_patients);

for (i in 1:n_patients) {
  int curr_patient_visit_pos, curr_patient_visit_end;
  (curr_patient_visit_pos, curr_patient_visit_end) = get_pos(patient_visit_pos, i);
  
  int first_visit = t_patient_visits[curr_patient_visit_pos];
  
  for (v in curr_patient_visit_pos:curr_patient_visit_end) {
    t_patient_visit_idx[v] = t_patient_visits[v] - first_visit + 1;
    
    if (t_patient_visits[v] <= 0) {
      n_patient_screening_visits[i] += 1;
    }
  }
  
  if (n_patient_screening_visits[i] == 0) {
    fatal_error("Patient ", i, " has no pre-screening visits.");
  }
}

// Time grid for full states computation: [1, 2, 3, ..., max_t_width]
// Each value represents weeks since patient's first visit
row_vector[max_t_width] time_since_first_visit = linspaced_row_vector(max_t_width, 1, max_t_width);

// Upper triangular matrix for cumulative sum integration of time-varying rates
// Column t contains sum of rates from columns 1 to t-1
// Used when enable_patient_process_noise_tr = 1
matrix[max_t_width, max_t_width] cumsum_integration_matrix = rep_matrix(0, max_t_width, max_t_width);
for (t in 2:max_t_width) {
  cumsum_integration_matrix[1:(t-1), t] = rep_vector(1, t-1);
}

// Combined flag: any process noise enabled (pop-level or patient-level)
int enable_any_process_noise_tr = enable_pop_process_noise_tr || enable_patient_process_noise_tr;

// ============================================================================
// FORECAST VISIT INFRASTRUCTURE
// ============================================================================
// These variables define the visit structure for forecasting and are shared
// across all biomarker modules (tumor, PSA, etc.)

int<lower = 1> last_predict_visit = max_all_t;

// Compute unique visits at population level (needed by state_space for visit_cumsum_mat)
int n_pop_unique_visits = num_unique(t_patient_visits, 0);
array[n_pop_unique_visits] int pop_unique_visits = unique(t_patient_visits, 0);

// Map patient visits to population unique visits
array[sum(n_patient_visits)] int<lower = 1> patient2pop_unique_visit_idx =
  get_level2level_idx(pop_unique_visits, t_patient_visits, patient_visit_pos);

// Last observed visit and forecast counts per patient
array[n_patients] int<lower = 1> patient_last_obs_visit = get_max_pos(t_patient_visits, patient_visit_pos);
array[n_patients] int<lower = 0, upper = last_predict_visit> n_patient_forecast_visits;

for (i in 1:n_patients) {
  n_patient_forecast_visits[i] = last_predict_visit - patient_last_obs_visit[i];
}

array[n_patients + 1] int<lower = 1> forecast_visits_pos = create_pos(n_patient_forecast_visits);

