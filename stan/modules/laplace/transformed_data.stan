// laplace/transformed_data.stan — Newton config and Laplace flat index arrays
//
// Routing index arrays (forecast_patient_idx, background_patient_idx, forecast_visit_pos,
// n_forecast_visits, n_background_patients) are declared in _base_transformed_data.stan
// because state_space/transformed_parameters.stan (shared across all models)
// uses them.

// ============================================================================
// LAPLACE NEWTON CONFIG + FLAT INDICES
// ============================================================================
// All laplace_n_* data variables are removed:
//   laplace_n_patients     -> n_background_patients (computed above)
//   laplace_n_total_visits -> not needed (patient_visit_pos[background_patient_idx[i]] used directly)
//   laplace_n_screening_visits[i] -> n_patient_screening_visits[background_patient_idx[i]]

real laplace_newton_tol = 1e-6;
int laplace_newton_max_iter = 4;  // Newton iterations for mode-finding; full Newton (mdivide_left_spd) gives quadratic convergence

// QR-projected covariates for non-target patients.
// Q_covar_design_matrix is now computed from all n_patients rows, so no R_full override needed.
matrix[n_background_patients, n_covar] laplace_Q_covar =
  enable_laplace_nontarget && n_covar > 0 && n_background_patients > 0
    ? Q_covar_design_matrix[background_patient_idx, :]
    : rep_matrix(0, n_background_patients, n_covar);

// Pre-compute flat indices for Laplace patients into the level intercept/slope arrays.
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_tr_intercept_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_frac_intercept_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_init_intercept_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_tr_slope_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_frac_slope_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_init_slope_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_ms_baseline_flat_idx;
array[n_background_patients, enable_laplace_nontarget ? n_levels : 0] int laplace_ms_slope_flat_idx;

if (enable_laplace_nontarget) {
  if (n_background_patients < 1)
    fatal_error("enable_laplace_nontarget=1 but n_background_patients < 1");

  for (i in 1:n_background_patients) {
    int p = background_patient_idx[i];  // unified patient index
    for (lv in 1:n_levels) {
      int group_id = patient_level_groups[p, lv];

      laplace_tr_intercept_flat_idx[i, lv] = enable_level_intercept_tr[lv]
        ? get_global_group_idx(enabled_level_pos_tr_intercept, lv, group_id) : 1;
      laplace_frac_intercept_flat_idx[i, lv] = enable_level_intercept_frac[lv]
        ? get_global_group_idx(enabled_level_pos_frac_intercept, lv, group_id) : 1;
      laplace_init_intercept_flat_idx[i, lv] = enable_level_intercept_init[lv]
        ? get_global_group_idx(enabled_level_pos_init_intercept, lv, group_id) : 1;

      laplace_tr_slope_flat_idx[i, lv] = enable_level_cov_tr[lv]
        ? get_global_group_idx(enabled_level_pos_tr_slope, lv, group_id) : 1;
      laplace_frac_slope_flat_idx[i, lv] = enable_level_cov_frac[lv]
        ? get_global_group_idx(enabled_level_pos_frac_slope, lv, group_id) : 1;
      laplace_init_slope_flat_idx[i, lv] = enable_level_cov_init[lv]
        ? get_global_group_idx(enabled_level_pos_init_slope, lv, group_id) : 1;

      laplace_ms_baseline_flat_idx[i, lv] = enable_ms_level_baseline_hazard[lv]
        ? get_global_group_idx(enabled_level_pos_ms_baseline, lv, group_id) : 1;
      laplace_ms_slope_flat_idx[i, lv] = enable_ms_level_cov[lv]
        ? get_global_group_idx(enabled_level_pos_ms_slope, lv, group_id) : 1;
    }
  }
}

// Index vehicle for reduce_sum: values 1..n_background_patients
// laplace_patient_range[i] = i, so the partial_sum can use start:end directly.
array[n_background_patients] int laplace_patient_range;
for (i in 1:n_background_patients)
  laplace_patient_range[i] = i;
