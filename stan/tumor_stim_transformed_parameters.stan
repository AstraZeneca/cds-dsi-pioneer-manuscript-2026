// Hazard ratio trial-level parameters 
vector[n_trials] tumor_stim_trial_intercept = add_trial_level ? raw_tumor_stim_trial_coef[, 1] * tumor_stim_trial_coef_sd[1] : rep_vector(0, n_trials);
matrix[n_trials, n_tumor_covar_col] tumor_stim_trial_coef;

for (s in 1:n_trials) {
  if (n_tumor_covar_col > 0) { 
    // TODO Add correlation between intercept and the coefs
    tumor_stim_trial_coef[s] = add_trial_level ? raw_tumor_stim_trial_coef[s, 2:] .* tumor_stim_trial_coef_sd[2:] : rep_row_vector(0, n_tumor_covar_col);
  }
}

// Hazard ratio organ-level parameters 
vector[n_tumor_locations] tumor_stim_location_intercept;
matrix[n_tumor_locations, n_tumor_covar_col] tumor_stim_location_coef;

tumor_stim_location_intercept = add_tumor_location_level ? raw_tumor_stim_location_coef[, 1] * tumor_stim_location_coef_sd[1] : rep_vector(0, n_tumor_locations);

if (n_tumor_covar_col > 0) { 
  for (l in 1:n_tumor_locations) {
    // TODO Add correlation between intercept and the coefs
    tumor_stim_location_coef[l] = 
      add_tumor_location_level ? raw_tumor_stim_location_coef[l, 2:] .* tumor_stim_location_coef_sd[2:] : rep_row_vector(0, n_tumor_covar_col);
  }
}
