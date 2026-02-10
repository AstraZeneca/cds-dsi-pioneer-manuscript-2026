// ============================================================================
// TUMOR-SPECIFIC FUNCTIONS
// ============================================================================
// These functions are specific to tumor/SLD measurements. They convert
// the general 2-component log-space state (from sf_state_space.stan) into
// tumor-specific quantities like Sum of Longest Diameters (SLD).
//
// For PSA or other biomarkers, create similar functions in modules/psa/functions.stan

/**
 * Calculate mean log(SLD) from two-component log-space states
 *
 * The two-component state represents:
 *   Component 1: Decreasing tumor burden (treatment effect)
 *   Component 2: Growing tumor burden (progression/resistance)
 *
 * Total SLD = exp(state[1]) + exp(state[2])
 * So log(SLD) = log_sum_exp(state[1], state[2])
 *
 * @param patient_states Matrix of states [n_visits × 2]
 * @param sum_tumor_size_baseline Baseline SLD measurement
 * @return Vector of log(SLD) means at each visit
 */
vector calc_log_sld_mean(matrix patient_states, real sum_tumor_size_baseline) {
  assert_equal(cols(patient_states), 2);

  return to_vector(log_sum_exp(patient_states[, 1], patient_states[, 2])) + log(sum_tumor_size_baseline);
}
