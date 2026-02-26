// ============================================================================
// TUMOR/SLD TRANSFORMED DATA
// ============================================================================
// SLD-specific preprocessing. Visit infrastructure (forecast_visits_pos, etc.)
// is now in _base_transformed_data.stan.

// Limit of detection for SLD measurements (cm)
real log_lod = log(0.1);

// RECIST response category constants
int CR = 1;
int PR = 2;
int SD = 3;
int PD = 4;

// Non-target lesion status constants
int NT_CR = 1;
int NT_PD = 2;

vector[sum(n_patient_visits)] log_sum_tumor_size = log(sum_tumor_size); // cm
vector[sum(n_patient_visits) - sum(n_patient_screening_visits)] post_treat_sld;

// Extract post-treatment SLD values
{
  int post_treat_pos = 1;

  for (i in 1:n_patients) {
    int curr_patient_visits_pos, curr_patient_visits_end;
    (curr_patient_visits_pos, curr_patient_visits_end) = get_pos(patient_visit_pos, i);

    for (m in curr_patient_visits_pos:curr_patient_visits_end) {
      if (t_patient_visits[m] > 0) {
        post_treat_sld[post_treat_pos] = sum_tumor_size[m];
        post_treat_pos += 1;
      }
    }
  }
}

// Array of measurement times used for GP modeling
array[max_all_t] real all_tumor_measure_t = linspaced_array(max_all_t, 1, max_all_t);

// ============================================================================
// BASELINE SLD AND NORMALIZATION
// ============================================================================
// Patient-level baseline SLD (in cm) for converting normalized states to absolute SLD
// Used by tumor likelihood and by other_events module for tumor covariates

vector[n_patients] log_baseline_sld;

// Normalize SLD by baseline for each patient (used by observation model)
vector<lower = 0>[sum(n_patient_visits)] normalized_sld;

{
  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);
    log_baseline_sld[i] = log(sum_tumor_size[visit_start]);
    normalized_sld[visit_start:visit_end] = sum_tumor_size[visit_start:visit_end] / sum_tumor_size[visit_start];
  }
}

// Biomarker-agnostic baseline for state_space module
vector[n_patients] baseline_obs_per_patient;
{
  for (i in 1:n_patients) {
    int visit_start = patient_visit_pos[i];
    baseline_obs_per_patient[i] = sum_tumor_size[visit_start];
  }
}
