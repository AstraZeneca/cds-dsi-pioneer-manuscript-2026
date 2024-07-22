// Per tumor population-level influence on hazard 
real<lower = 0> tumor_stim_pop_intercept; // DO NOT REMOVE; this is a per tumor intercept and not per patient intercept which is included in lambda.
row_vector[n_tumor_covar_col] tumor_stim_pop_coef;

// Trial level hierarchical tumor effect
matrix[add_trial_level ? n_trials : 0, n_tumor_covar_col + 1] raw_tumor_stim_trial_coef;
row_vector<lower = 0>[add_trial_level ? n_tumor_covar_col + 1 : 0] tumor_stim_trial_coef_sd;

// Organ level hierarchical tumor effect 
matrix[add_tumor_location_level ? n_tumor_locations : 0, n_tumor_covar_col + 1] raw_tumor_stim_location_coef;
row_vector<lower = 0>[add_tumor_location_level ? n_tumor_covar_col + 1 : 0] tumor_stim_location_coef_sd;  