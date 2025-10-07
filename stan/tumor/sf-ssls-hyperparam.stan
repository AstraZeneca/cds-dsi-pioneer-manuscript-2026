// GP parameters
real<lower = 0> pop_tumor_gp_rho_meanlog;
real<lower = 0> pop_tumor_gp_rho_sdlog;
real<lower = 0> log_patient_tumor_gp_rho_sd_sd;

// Process noise parameters
real<lower = 0> pop_decrease_process_sd_sd;
real<lower = 0> pop_growth_process_sd_sd;
real<lower = 0> process_corr_param;
real<lower = 0> measure_sd_sd;

// Population rate parameters (total + fraction)
real pop_log_total_rate_mean;
real<lower = 0> pop_log_total_rate_sd;
real pop_decrease_frac_logit_mean;
real<lower = 0> pop_decrease_frac_logit_sd;
real<lower = 0> trial_log_total_rate_sd_sd;
real<lower = 0> patient_log_total_rate_sd_sd;
real<lower = 0> patient_decrease_frac_logit_sd_sd; 

// (Removed: pop_log_net_rate_mean, pop_log_net_rate_sd,
//           pop_log_rate_ratio_mean, pop_log_rate_ratio_sd,
//           trial_log_net_rate_sd_sd, patient_log_net_rate_sd_sd,
//           patient_log_rate_ratio_sd_sd)

// Growth lag parameters
real growth_lag_mean;
real<lower = 0> growth_lag_sd;
real<lower = 0> patient_log_growth_lag_sd_sd;
real<lower = 0> log_growth_transition_rate_sd;

// Correlation parameters
real<lower = 0> rate_corr_param;

// Proportion parameters
real pop_decrease_prop_logis_mean;
real<lower = 0> pop_decrease_prop_logis_sd;
real<lower = 0> trial_decrease_prop_logis_sd_sd;
real<lower = 0> patient_decrease_prop_logis_sd_sd;

real<lower = 0> log_lod_sd;

int<lower = 1> n_causes; // Death and non-target PD

vector<lower = 0>[n_causes] log_lambda_gp_pop_alpha_sd;
vector<lower = 0>[n_causes] log_lambda_gp_pop_rho_alpha, log_lambda_gp_pop_rho_beta;
vector[n_causes] log_lambda_gp_pop_intercept_mean;
vector<lower = 0>[n_causes] log_lambda_gp_pop_intercept_sd;

vector<lower = 0>[n_causes] log_lambda_gp_trial_alpha_sd;
vector<lower = 0>[n_causes] log_lambda_gp_trial_rho_alpha, log_lambda_gp_trial_rho_beta;
vector<lower = 0>[n_causes] log_lambda_gp_trial_intercept_sd_sd;

// Prior hyperparameters for covariate effects
// Shift linear model from total rate to fraction (decrease share) on logit scale
vector[n_covar] pop_decrease_frac_logit_coef_mean;
vector<lower=0>[n_covar] pop_decrease_frac_logit_coef_sd;
row_vector<lower=0>[n_covar] trial_decrease_frac_logit_coef_sd_sd;
// Patient-level (original beta scale) per-covariate SD hyperparameters for decrease fraction logit coefficients
row_vector<lower=0>[n_covar] patient_decrease_frac_logit_coef_sd_sd;

// Initial state proportion covariate effects remain
vector[n_covar] pop_decrease_prop_logis_coef_mean;
vector<lower=0>[n_covar] pop_decrease_prop_logis_coef_sd;