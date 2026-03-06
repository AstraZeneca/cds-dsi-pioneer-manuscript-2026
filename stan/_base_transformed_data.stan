#include "_base_hierarchy_transformed_data.stan"

array[n_patients + 1] int<lower = 1, upper = sum(n_patient_visits) + 1> patient_visit_pos = create_pos(n_patient_visits);

int max_all_t = max(max(t_patient_visits) + 1, extend_max_all_t); // Latest measurement time or extended time, whichever is greater
int<lower = 0> max_t_width = max_all_t - min(t_patient_visits) + 1;


array[sum(n_patient_visits)] int<lower = 1> t_patient_visit_idx; // Index of each patient visit relative to the first visit for each patient 

array[n_patients] int<lower = 0> n_patient_screening_visits = zeros_int_array(n_patients);

for (i in 1:n_patients) {
  int curr_patient_visit_pos, curr_patient_visit_end;
  (curr_patient_visit_pos, curr_patient_visit_end) = get_pos(patient_visit_pos, i);
  
  // Count screening visits first to locate the baseline.
  for (v in curr_patient_visit_pos:curr_patient_visit_end) {
    if (t_patient_visits[v] <= 0) {
      n_patient_screening_visits[i] += 1;
    }
  }

  if (n_patient_screening_visits[i] == 0) {
    fatal_error("Patient ", i, " has no pre-screening visits.");
  }

  // Baseline = last screening visit (latest visit with week <= 0).
  int baseline_visit_idx = curr_patient_visit_pos + n_patient_screening_visits[i] - 1;
  int baseline_week = t_patient_visits[baseline_visit_idx];

  // Anchor t_patient_visit_idx at the baseline so that t=1 corresponds to
  // the baseline for all patients. Pre-screening visits get idx <= 0 but are
  // clamped to 1 since they are never used in the likelihood.
  for (v in curr_patient_visit_pos:curr_patient_visit_end) {
    int raw_idx = t_patient_visits[v] - baseline_week + 1;
    t_patient_visit_idx[v] = raw_idx > 0 ? raw_idx : 1;
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

