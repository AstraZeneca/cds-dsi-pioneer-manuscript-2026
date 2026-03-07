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
