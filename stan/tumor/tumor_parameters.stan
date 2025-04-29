real<lower = 0> pop_tumor_intercept;

// Multilevel intercepts

real<lower = 0> trial_tumor_intercept_sd;
vector[add_trial_level_tumor_intercept ? n_trials : 0] raw_trial_tumor_intercept_effect;

real<lower = 0> patient_tumor_intercept_sd;
vector[n_patients] raw_patient_tumor_intercept_effect;

real<lower = 0> tumor_measure_error_sd; 
// real<lower = 0> pop_lod;