// ============================================================================
// Visit Infrastructure Transformed Data
// Shared by all models (full biomarker models and standalone multistate).
//
// Requires in scope:
//   n_patients         — from _hierarchy_data.stan
//   n_patient_visits[] — from _visit_data.stan
//   t_patient_visits[] — from _visit_data.stan
//   max_all_t          — from data (standalone) or derived inline (full models)
// ============================================================================

array[n_patients + 1] int<lower=1, upper=sum(n_patient_visits) + 1> patient_visit_pos = create_pos(n_patient_visits);
int n_total_visits = sum(n_patient_visits);

// Count screening visits per patient (visits with t <= 0)
array[n_patients] int<lower=0> n_patient_screening_visits = zeros_int_array(n_patients);
for (i in 1:n_patients) {
  int v_start, v_end;
  (v_start, v_end) = get_pos(patient_visit_pos, i);
  for (v in v_start:v_end) {
    if (t_patient_visits[v] <= 0) n_patient_screening_visits[i] += 1;
  }
  if (n_patient_screening_visits[i] == 0) {
    fatal_error("Patient ", i, " has no pre-screening visits.");
  }
}

// Forecast visit infrastructure
int<lower=1> last_predict_visit = max_all_t;
array[n_patients] int<lower=1> patient_last_obs_visit = get_max_pos(t_patient_visits, patient_visit_pos);
array[n_patients] int<lower=0> n_patient_forecast_visits;
for (i in 1:n_patients) {
  n_patient_forecast_visits[i] = last_predict_visit - patient_last_obs_visit[i];
}
array[n_patients + 1] int<lower=1> forecast_visits_pos = create_pos(n_patient_forecast_visits);
