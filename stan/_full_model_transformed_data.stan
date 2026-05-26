// ============================================================================
// Full Biomarker Model Transformed Data
//
// Requires in scope from preceding includes/inline derivations:
//   patient_visit_pos[], n_patient_screening_visits[]  — from _visit_transformed_data.stan
//   max_all_t, max_t_width                             — derived inline in model file
//   n_forecast_patients, forecast_split_level, forecast_group — from _full_model_data.stan
//   n_patients, patient_trial, patient_level_groups    — from _hierarchy_data.stan
//   trial_patient_pos                                  — from _hierarchy_transformed_data.stan
// ============================================================================

// Index of each patient visit relative to the first visit for each patient
array[sum(n_patient_visits)] int<lower = 1> t_patient_visit_idx;

for (i in 1:n_patients) {
  int curr_patient_visit_pos, curr_patient_visit_end;
  (curr_patient_visit_pos, curr_patient_visit_end) = get_pos(patient_visit_pos, i);

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

// Combined flag: any process noise enabled (pop-level or patient-level),
// and whether the full states grid [n_patients × max_t_width] is needed:
//   - process noise: rates vary by timepoint, need states at every week
//   - enable_states_full_grid: explicit flag (future-proofing)
//   - ungated continuous time-varying covariate for 0->1 multistate transition
//   - gated+latent 0->1, 0->2, or 0->3 TV cov (requires grid when process noise ON)
// ms_needs_inline_burden: when process noise OFF but MS needs modeled burden,
//   compute burden trajectory analytically per-patient instead of building full grid
int enable_any_process_noise_tr;
int need_states_full_grid;
int ms_needs_inline_burden;
(enable_any_process_noise_tr, need_states_full_grid, ms_needs_inline_burden) = compute_full_model_grid_flags(
  enable_pop_process_noise_tr, enable_patient_process_noise_tr,
  enable_states_full_grid,
  enable_ms_pop_time_varying_cov, n_time_varying_covar,
  enable_ms_01, enable_ms_visit_gated_01, enable_ms_visit_gated_latent_01,
  enable_ms_02_time_varying_cov, enable_ms_03_time_varying_cov
);
print("need_states_full_grid = ", need_states_full_grid,
      " (process_noise=", enable_any_process_noise_tr,
      ", ms_needs_inline_burden=", ms_needs_inline_burden, ")");

// ============================================================================
// POPULATION UNIQUE VISITS
// ============================================================================

// Compute unique visits at population level (needed by state_space for visit_cumsum_mat)
int n_pop_unique_visits = num_unique(t_patient_visits, 0);
array[n_pop_unique_visits] int pop_unique_visits = unique(t_patient_visits, 0);

// Map patient visits to population unique visits
array[sum(n_patient_visits)] int<lower = 1> patient2pop_unique_visit_idx =
  get_level2level_idx(pop_unique_visits, t_patient_visits, patient_visit_pos);

// ============================================================================
// PATIENT ROUTING: Forecast vs Background index arrays
// ============================================================================
// Derived from n_forecast_patients, forecast_split_level, forecast_group
// (all declared in _full_model_data.stan).
// Lives here (not laplace/transformed_data.stan) because state_space and
// multistate transformed_parameters are shared across all models.
//
// Forecast patients = trial patients we generate quantities for
// Background patients = RWD/non-target (contribute to likelihood, no GQ)
//
// When forecast_split_level = 0:
//   forecast_patient_idx = [1, 2, ..., n_patients]  (all patients are forecast)
//   background_patient_idx is empty

int n_background_patients = n_patients - n_forecast_patients;
array[n_forecast_patients] int forecast_patient_idx;
array[n_background_patients] int background_patient_idx;

{
  int j_fc = 0;
  int j_bg = 0;
  for (p in 1:n_patients) {
    if (forecast_split_level > 0 &&
        patient_level_groups[p, forecast_split_level] != forecast_group) {
      j_bg += 1;
      background_patient_idx[j_bg] = p;
    } else {
      j_fc += 1;
      forecast_patient_idx[j_fc] = p;
    }
  }
  if (j_fc != n_forecast_patients)
    fatal_error("n_forecast_patients=", n_forecast_patients,
                " but found ", j_fc, " forecast patients via routing key");
}

// Forecast-local visit position array — parallel to patient_visit_pos but for forecast patients only.
int n_forecast_visits = sum(n_patient_visits[forecast_patient_idx]);
array[n_forecast_patients + 1] int forecast_visit_pos = create_pos(n_patient_visits[forecast_patient_idx]);

// Forecast-patient position array: like trial_patient_pos but counts only forecast patients
// per trial. Used to slice n_forecast_patients-sized endpoint arrays in GQ by trial.
array[n_trials + 1] int forecast_trial_patient_pos;
{
  array[n_trials] int n_forecast_per_trial = rep_array(0, n_trials);
  for (j in 1:n_forecast_patients) {
    n_forecast_per_trial[patient_trial[forecast_patient_idx[j]]] += 1;
  }
  forecast_trial_patient_pos = create_pos(n_forecast_per_trial);
}

if (forecast_visit_pos[n_forecast_patients + 1] - 1 != n_forecast_visits)
  fatal_error("n_forecast_visits=", n_forecast_visits, " inconsistent with forecast_visit_pos sum");

// ============================================================================
// DUAL INDEXING CONVENTION
// ============================================================================
// Two visit position systems coexist after this point:
//
//   UNIFIED positions  — patient_visit_pos[p]  where p = forecast_patient_idx[j]
//                                                    or p = background_patient_idx[i]
//     Used for: all data arrays (normalized_psa, t_patient_visit_idx,
//               ms_final_state, ms_time_*, psa_values, log_baseline_psa, ...)
//
//   FORECAST-LOCAL positions — forecast_visit_pos[j]  where j = 1..n_forecast_patients
//     Used for: states[n_forecast_visits, 2] only
//
// Pattern in every forecast patient loop:
//   for (j in 1:n_forecast_patients) {
//     int p = forecast_patient_idx[j];            // unified patient index
//     int data_start = patient_visit_pos[p];      // → data arrays
//     int state_start = forecast_visit_pos[j];    // → states matrix
//   }
