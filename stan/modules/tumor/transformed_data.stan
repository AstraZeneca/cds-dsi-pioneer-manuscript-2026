// ============================================================================
// TUMOR/SLD TRANSFORMED DATA
// ============================================================================
// SLD-specific preprocessing. Visit infrastructure (forecast_visits_pos, etc.)
// is now in _base_transformed_data.stan.

// Limit of detection for SLD measurements (cm)
real log_lod = log(0.1);

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

// Array of time points used for GP modeling in the multistate module
array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

// ============================================================================
// BASELINE SLD AND NORMALIZATION
// ============================================================================
// Patient-level baseline SLD (in cm) for converting normalized states to absolute SLD
// Used by tumor likelihood and by multistate module for tumor covariates

vector[n_patients] log_baseline_sld;

// ============================================================================
// SLD NORMALIZATION CONSTANTS FOR PROPORTIONAL HAZARDS
// ============================================================================
// Robust normalization statistics (median and IQR) from observed log(SLD)
// Using median/IQR instead of mean/SD for robustness to outliers
// Used by multistate module for time-varying covariates

real median_log_sld_obs;
real iqr_log_sld_obs;

// Compute robust normalization statistics from observed log(SLD)
{
  // Get all observed log(SLD) values across all patients and visits
  // Filter out zeros to avoid -Inf
  vector[sum(n_patient_visits)] log_sld_all_obs;
  int n_positive = 0;

  for (i in 1:sum(n_patient_visits)) {
    if (sum_tumor_size[i] > 0) {
      n_positive += 1;
      log_sld_all_obs[n_positive] = log(sum_tumor_size[i]);
    }
  }

  // Compute robust normalization using median and IQR from positive observations
  {
    array[3] real quantiles_obs = quantile(log_sld_all_obs[1:n_positive], {0.25, 0.5, 0.75});
    median_log_sld_obs = quantiles_obs[2];
    iqr_log_sld_obs = quantiles_obs[3] - quantiles_obs[1];
  }
}

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
