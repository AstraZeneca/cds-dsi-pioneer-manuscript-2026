vector<lower = 0>[n_tumor_separate_trials] tumor_mean; // Actual tumor mean, not the lognormal mean
vector<lower = 0>[n_tumor_separate_trials] tumor_sd; // Actual tumor SD, not the lognormal one. Homoskedastic for now.
// real<lower = 0> pop_lod;

vector<lower = 0>[n_tumor_separate_trials] pop_tumor_gp_alpha;
vector<lower = 0>[n_tumor_separate_trials] pop_tumor_gp_rho;
array[n_tumor_separate_trials] vector[model_all_measures ? max_t_width : n_pop_unique_visits] latent_pop_tumor_gp_eta;


// Multilevel intercepts

// real<lower = 0> trial_tumor_gp_intercept_sd;
// vector[add_trial_level_tumor_gp ? n_trials : 0] trial_tumor_gp_intercept_effect;
// 
// real<lower = 0> patient_tumor_gp_intercept_sd;
// vector[add_patient_level_tumor_gp ? n_patients : 0] patient_tumor_gp_intercept_effect;

// vector<lower = 0>[use_tumor_model && multilevel_tumor ? n_patients : 0] tumor_gp_intercept_sd;
// vector[use_tumor_model && multilevel_tumor ? sum(n_patient_tumors) : 0] tumor_gp_intercept_effect;
