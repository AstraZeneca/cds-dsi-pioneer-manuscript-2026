int<lower = 0, upper = 1> fit_tumor_data; // Sample of prior prediction only
int<lower = 0, upper = 1> predict_missing_sizes; 
int<lower = 0, upper = 1> gen_tumor_sizes; // Generate data given prior/posterior of model parameters
int<lower = 0, upper = 1> model_all_measures; // Doesn't do anything right now

int<lower = 0, upper = 1> add_trial_level_tumor_intercept;
int<lower = 0, upper = 1> add_trial_level_tumor_gp_param;

// Hyperparam
real tumor_mean_mean;
real<lower = 0> tumor_mean_sd;
real<lower = 0> tumor_measure_error_sd_sd;

real<lower = 0> patient_tumor_intercept_sd_sd;
real<lower = 0> trial_tumor_intercept_sd_sd;

real<lower = 0> pop_tumor_gp_alpha_sd;
real<lower = 0> pop_tumor_gp_rho_meanlog;
real<lower = 0> pop_tumor_gp_rho_sdlog;

real<lower = 0> log_trial_tumor_gp_rho_sd_sd;
real<lower = 0> log_patient_tumor_gp_rho_sd_sd;