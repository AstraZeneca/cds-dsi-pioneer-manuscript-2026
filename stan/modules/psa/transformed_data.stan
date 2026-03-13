// ============================================================================
// PSA NORMALIZATION AND PREPROCESSING
// ============================================================================

vector[sum(n_patient_visits)] log_psa_values;
array[n_patients] real log_baseline_psa;

// Normalized PSA (relative to baseline) - only for measured values
vector[sum(n_patient_visits)] normalized_psa;

// Extract post-treatment PSA values (parallel to post_treat_sld in tumor module)
vector[sum(n_patient_visits) - sum(n_patient_screening_visits)] post_treat_psa;

// Compute normalized PSA per patient
{
  int post_treat_pos = 1;

  for (i in 1:n_patients) {
    int visit_start, visit_end;
    (visit_start, visit_end) = get_pos(patient_visit_pos, i);

    // Baseline is the LAST screening visit (week 0, ADaM ABLFL=Y).
    // Treatment visits (week > 0) start at visit_start + n_screening.
    int baseline_idx = visit_start + n_patient_screening_visits[i] - 1;
    real baseline_psa = psa_values[baseline_idx];
    log_baseline_psa[i] = log(baseline_psa);


    for (v in visit_start:visit_end) {
      if (psa_measured[v]) {
        normalized_psa[v] = psa_values[v] / baseline_psa;
        log_psa_values[v] = log(psa_values[v]);
      } else {
        normalized_psa[v] = 0;  // Placeholder for unmeasured
        log_psa_values[v] = 0;
      }

      // Extract post-treatment values
      if (t_patient_visits[v] > 0) {
        post_treat_psa[post_treat_pos] = psa_values[v];
        post_treat_pos += 1;
      }
    }
  }
}

// Biomarker-agnostic baseline for state_space module
vector[n_patients] baseline_obs_per_patient;
{
  for (i in 1:n_patients) {
    baseline_obs_per_patient[i] = exp(log_baseline_psa[i]);
  }
}

// ============================================================================
// GP TIME GRID
// ============================================================================
// Array of absolute time points used for GP modeling in the multistate module.
array[max_all_t] real all_measure_t = linspaced_array(max_all_t, 1, max_all_t);

// ============================================================================
// PSA NORMALIZATION CONSTANTS FOR PROPORTIONAL HAZARDS
// ============================================================================
// Robust normalization statistics (median and IQR) from observed log(PSA)
// Used by multistate module for time-varying covariates
// Only computed from measured PSA values (psa_measured == 1 and non-zero)

real median_log_psa_obs;
real iqr_log_psa_obs;

{
  int n_positive = 0;
  for (i in 1:sum(n_patient_visits)) {
    if (psa_measured[i] && log_psa_values[i] != 0) n_positive += 1;
  }
  vector[n_positive] obs_log_psa;
  int obs_idx = 1;
  for (i in 1:sum(n_patient_visits)) {
    if (psa_measured[i] && log_psa_values[i] != 0) {
      obs_log_psa[obs_idx] = log_psa_values[i];
      obs_idx += 1;
    }
  }
  array[3] real quantiles_obs = quantile(obs_log_psa, {0.25, 0.5, 0.75});
  median_log_psa_obs = quantiles_obs[2];
  iqr_log_psa_obs = quantiles_obs[3] - quantiles_obs[1];
}
