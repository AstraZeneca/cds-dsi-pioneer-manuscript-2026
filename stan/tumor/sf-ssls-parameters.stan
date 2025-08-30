real pop_log_net_rate;            // Population-level net rate (log(d-g))
real<lower = 0> pop_log_rate_ratio; // Population-level ratio (log(d/g))

// Trial-level variation for net rate only
real<lower=0> trial_log_net_rate_sd;
vector[pop_rates_param_only || !add_trial_level_net_rate ? 0 : n_train_trials] raw_trial_log_net_rate;

// Patient-level variation for net rate only
real<lower=0> patient_log_net_rate_sd;
// vector<offset = pop_log_net_rate, multiplier = patient_log_net_rate_sd>[n_patients] patient_log_net_rate;
vector[pop_rates_param_only ? 0 : n_train_patients] raw_patient_log_net_rate;

// real<lower=0> patient_log_rate_ratio_sd;
// vector[pop_rates_param_only ? 0 : n_train_patients] raw_patient_log_rate_ratio; 

real pop_log_growth_lag;
real pop_log_growth_transition_rate;

real<lower = 0> patient_log_growth_lag_sd;
// vector<offset = pop_log_growth_lag, multiplier = patient_log_growth_lag_sd>[n_patients] patient_log_growth_lag;
vector[pop_growth_lag_param_only ? 0 : n_train_patients] raw_patient_log_growth_lag;

// Patient-level GP
// row_vector<lower = 0>[2] pop_tumor_gp_alpha;
real log_pop_tumor_gp_rho;

// Hierarchical length-scale parameter
// vector<lower = 0>[2] log_trial_tumor_gp_rho_sd;
// vector[add_trial_level_tumor_gp_param ? n_trials : 0] raw_log_trial_tumor_gp_rho_effect; 

real<lower = 0> log_patient_tumor_gp_rho_sd;
vector[independ_long_process_noise || pop_rho_param_only ? 0 : n_train_patients] raw_log_patient_tumor_gp_rho_effect; 

matrix[n_total_train_visits_m1, 2] raw_patient_process_noise;

// real<lower=0> pop_decrease_process_sd;
// real<lower=0> pop_growth_process_sd;
vector<lower = 0>[2] pop_process_sd;
cholesky_factor_corr[independ_cross_process_noise ? 0 : 2] L_process_corr;

real<lower=0> measure_sd;        // Measurement noise variance

// real<lower = 0> patient_log_decrease_rate_sd;
// real<lower = 0> patient_log_growth_rate_sd;

real pop_decrease_prop_logis;

real<lower = 0> trial_decrease_prop_logis_sd;
vector[pop_initial_states_param_only || !add_trial_level_prop ? 0 : n_train_trials] raw_trial_decrease_prop_logis;

real<lower = 0> patient_decrease_prop_logis_sd;
// vector<offset = pop_decrease_prop_logis, multiplier = patient_decrease_prop_logis_sd>[n_patients] patient_decrease_prop_logis;
vector[pop_initial_states_param_only ? 0 : n_train_patients] raw_patient_decrease_prop_logis;

// real log_lod;

// Covariate effects on rates
vector[n_covar] pop_log_net_rate_coef;      // Population-level covariate effects on net rate (original space)
vector[n_covar] pop_decrease_prop_logis_coef;

// Optional: hierarchical covariate effects
row_vector<lower=0>[pop_covar_coef_only ? 0 : n_covar] trial_log_net_rate_coef_sd;
matrix[pop_covar_coef_only ? 0 : n_train_trials, n_covar] raw_trial_log_net_rate_coef;

// real<lower=0> patient_log_net_rate_coef_sd;
// matrix[pop_rates_param_only ? 0 : n_train_patients, n_covar] raw_patient_log_net_rate_coef;
