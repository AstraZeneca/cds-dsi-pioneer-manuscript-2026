vector<lower = 0>[n_tumor_separate_trials] tumor_mean;
vector<lower = 0>[n_tumor_separate_trials] tumor_measure_error_sd; 
// real<lower = 0> pop_lod;

// Multilevel intercepts

real<lower = 0> trial_tumor_gp_intercept_sd;
vector[add_trial_level_tumor ? n_trials : 0] raw_trial_tumor_gp_intercept_effect;

vector<lower = 0>[n_tumor_separate_trials] patient_tumor_gp_intercept_sd;
vector[n_patients] raw_patient_tumor_gp_intercept_effect;

// GP
// Population-level additative GP
vector<lower = 0>[patient_gp_only ? 0 : n_tumor_separate_trials] pop_tumor_gp_alpha;
vector<lower = 0>[patient_gp_only ? 0 : n_tumor_separate_trials] pop_tumor_gp_rho;
vector[patient_gp_only ? 0 : (separate_trial_tumor_gp ? sum(n_trial_unique_visits) : n_pop_unique_visits)] pop_tumor_gp_eta;

// Trial-level additative GP 
real<lower = 0> trial_tumor_gp_alpha;
real<lower = 0> trial_tumor_gp_rho;
vector[add_trial_level_tumor && !patient_gp_only ? sum(n_trial_unique_visits) : 0] trial_tumor_gp_eta;

// Patient-level GP
vector<lower = 0>[n_tumor_separate_trials] patient_tumor_gp_alpha;
// If we're using a multilevel model of rho we're actually in the log space.
vector<lower = (multilevel_gp_param ? negative_infinity() : 0)>[n_tumor_separate_trials] patient_tumor_gp_rho;
vector[sum(n_patient_unique_visits)] patient_tumor_gp_eta;

// Hierarchical length-scale parameter
real<lower = 0> log_trial_rho_sd;
vector[multilevel_gp_param ? n_trials : 0] raw_log_trial_rho; 

real<lower = 0> log_patient_rho_sd;
vector[multilevel_gp_param ? n_patients : 0] raw_log_patient_rho; 
