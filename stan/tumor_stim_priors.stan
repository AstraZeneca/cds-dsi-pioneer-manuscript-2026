tumor_stim_pop_intercept ~ normal(0, tumor_stim_pop_intercept_sd);

if (add_trial_level) { 
  tumor_stim_trial_coef_sd ~ normal(0, tumor_stim_trial_coef_sd_sd[:(n_tumor_covar_col + 1)]);
  
  to_vector(raw_tumor_stim_trial_coef) ~ std_normal();
}

if (add_tumor_location_level) { 
  tumor_stim_location_coef_sd ~ normal(0, tumor_stim_location_coef_sd_sd[:(n_tumor_covar_col + 1)]);
  to_vector(raw_tumor_stim_location_coef) ~ std_normal();
}

if ((tumor_hazard_type > 0 && tumor_hazard_type < 4) || tumor_hazard_type > 4) {
  tumor_stim_pop_coef[1] ~ normal(0, tumor_stim_pop_coef_sd[1]);

  if (tumor_hazard_type == 1 || tumor_hazard_type == 2) {
    tumor_stim_pop_coef[2] ~ normal(0, tumor_stim_pop_coef_sd[2]);

    if (tumor_hazard_type == 2 || tumor_hazard_type == 5) {
      tumor_stim_pop_coef[3] ~ normal(0, tumor_stim_pop_coef_sd[3]);
    }
    
    if (tumor_hazard_type == 5) {
      tumor_stim_pop_coef[4] ~ normal(0, tumor_stim_pop_coef_sd[4]);
      tumor_stim_pop_coef[5] ~ normal(0, tumor_stim_pop_coef_sd[5]);
    }
  }
}
