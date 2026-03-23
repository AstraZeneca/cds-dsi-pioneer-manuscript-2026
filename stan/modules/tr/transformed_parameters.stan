// tr/transformed_parameters.stan — active linear predictor assembly for total rate module
// Optimized: uses pre-computed flat indices from transformed_data for direct gather

// Population covariate effects
vector[n_hmc_patients] tr_linpred_pop = enable_pop_cov_tr ?
  (Q_covar_design_matrix[hmc_patient_idx, :] * tr_coef_qr_pop) : zeros_vector(n_hmc_patients);

// ===== INTERCEPT EFFECTS =====
// Step 1: Scale all raw effects at once (vectorized per level)
vector[n_enabled_groups_tr_intercept] tr_scaled_level_intercept;
for (lv in 1:n_levels) {
  if (enable_level_intercept_tr[lv]) {
    int lv_start, lv_end;
    (lv_start, lv_end) = get_pos(enabled_level_pos_tr_intercept, lv);
    tr_scaled_level_intercept[lv_start:lv_end] =
      tr_sd_level_intercept[lv] * tr_raw_level_intercept[lv_start:lv_end];
  }
}

// Step 2: Gather using pre-computed flat indices (no intermediate array creation)
vector[n_hmc_patients] tr_linpred_level_intercepts = zeros_vector(n_hmc_patients);
for (lv in 1:n_levels) {
  if (enable_level_intercept_tr[lv]) {
    tr_linpred_level_intercepts += tr_scaled_level_intercept[patient_tr_intercept_flat_idx[hmc_patient_idx, lv]];
  }
}

// ===== COVARIATE SLOPE EFFECTS =====
// Step 1: Scale all raw slope effects at once
matrix[n_enabled_groups_tr_slope, n_covar] tr_scaled_level_slope;
if (n_covar > 0 && n_enabled_groups_tr_slope > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_tr[lv]) {
      int lv_start, lv_end;
      (lv_start, lv_end) = get_pos(enabled_level_pos_tr_slope, lv);
      tr_scaled_level_slope[lv_start:lv_end, :] =
        tr_raw_level_slope[lv_start:lv_end, :] .*
        rep_matrix(tr_sd_level_slope[lv]', lv_end - lv_start + 1);
    }
  }
}

// Step 2: Gather and compute dot products using pre-computed flat indices
vector[n_hmc_patients] tr_linpred_level_slopes = zeros_vector(n_hmc_patients);
if (n_covar > 0) {
  for (lv in 1:n_levels) {
    if (enable_level_cov_tr[lv]) {
      tr_linpred_level_slopes += rows_dot_product(
        Q_covar_design_matrix[hmc_patient_idx, :],
        tr_scaled_level_slope[patient_tr_slope_flat_idx[hmc_patient_idx, lv], :]
      );
    }
  }
}

// ===== FINAL LINEAR PREDICTOR =====
vector[n_hmc_patients] tr_loc_patient = tr_loc_pop
  + tr_linpred_pop
  + tr_linpred_level_intercepts
  + tr_linpred_level_slopes; 

vector[enable_patient_process_noise_tr ? n_hmc_patients : 0] tr_log_sd_patient_process_noise;
vector[enable_patient_process_noise_tr ? n_hmc_patients : 0] tr_phi_patient_process_noise;
matrix[enable_patient_process_noise_tr ? n_hmc_patients : 0, max_t_width] tr_patient_process_noise;

if (enable_patient_process_noise_tr) {
  // Non-centered parameterization: work in log-space, exponentiate once at the end
  // If patient hierarchy is disabled, all patients get population value
  if (enable_patient_process_noise_sd_tr) {
    tr_log_sd_patient_process_noise = tr_log_sd_pop_process_noise[1] + tr_sd_patient_log_sd_process_noise[1] * tr_raw_patient_log_sd_process_noise;
  } else {
    tr_log_sd_patient_process_noise = rep_vector(tr_log_sd_pop_process_noise[1], n_hmc_patients);
  }

  // Map to [0,1] via inv_logit link function
  // Population phi is inv_logit(logit_phi_pop), patient deviations on logit scale
  if (enable_patient_process_noise_phi_tr) {
    tr_phi_patient_process_noise = inv_logit(
      tr_logit_phi_pop_process_noise[1] + tr_sd_patient_phi_process_noise[1] * tr_raw_patient_phi_process_noise
    );
  } else {
    tr_phi_patient_process_noise = rep_vector(inv_logit(tr_logit_phi_pop_process_noise[1]), n_hmc_patients);
  }

  vector[n_hmc_patients] tr_sd_patient_process_noise = exp(tr_log_sd_patient_process_noise);

  // AR(1) process: deviation[t] = phi * deviation[t-1] + sigma * innovation[t]
  // This creates time-varying deviations with mean 0 that revert to baseline
  tr_patient_process_noise[, 1] = tr_raw_patient_process_noise[, 1] .* tr_sd_patient_process_noise;

  for (t in 2:max_t_width) {
    tr_patient_process_noise[, t] = tr_phi_patient_process_noise .* tr_patient_process_noise[, t - 1]
      + tr_raw_patient_process_noise[, t] .* tr_sd_patient_process_noise;
  }
}

// Population-level time-varying process noise (shared AR(1) across all patients)
row_vector[enable_pop_process_noise_tr ? max_t_width : 0] tr_pop_process_noise;

if (enable_pop_process_noise_tr) {
  real sigma_pop = exp(tr_log_sd_pop_process_noise_pop[1]);
  // log(inv_logit(x)) = x - log1p_exp(x), more stable than log(inv_logit(x))
  real log_phi_pop = tr_logit_phi_pop_process_noise_pop[1] - log1p_exp(tr_logit_phi_pop_process_noise_pop[1]);

  // Compute φ^t for t = 1..max_t_width (vectorized via exp/log)
  row_vector[max_t_width] phi_powers = exp(log_phi_pop * linspaced_row_vector(max_t_width, 1, max_t_width));

  // Scale raw innovations
  row_vector[max_t_width] z_pop = sigma_pop * tr_raw_pop_process_noise;

  // Vectorized AR(1): y[t] = Σ_{k=1}^{t} φ^(t-k) * z[k] = φ^t * cumsum(z/φ^k)
  tr_pop_process_noise = cumulative_sum(z_pop ./ phi_powers) .* phi_powers;
}