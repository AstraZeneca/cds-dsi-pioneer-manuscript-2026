// Hyperparam for baseline hazards
array[separate_baseline_hazard ? n_trials : 1] vector[2] log_crcr_lambda_gp_intercept_mean;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_crcr_lambda_gp_intercept_sd;

// Multilevel SD for baseline hazard intercept
real<lower = 0> log_crcr_lambda_gp_trial_intercept_sd_sd;

vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_crcr_lambda_gp_alpha_sd;

// Multilevel SD for alpha only. Reusing the the rho hyperparam for the multilevel GP.
real<lower = 0> log_crcr_lambda_gp_trial_alpha_sd;

// Hyperparam for GP rho 
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_crcr_lambda_gp_rho_alpha;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_crcr_lambda_gp_rho_beta;

// Hyperparam for proportional hazard parameters
array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[n_tumor_covar] crcr_tumor_stim_pop_coef_sd;
array[separate_prop_hazard ? n_trials : 1] vector[n_covar] crcr_covar_effect_mean;
array[separate_prop_hazard ? n_trials : 1] vector<lower = 0>[n_covar] crcr_covar_effect_sd;

// Multilevel param for prop hazard parameters
real<lower = 0> crcr_covar_trial_sd_sd;
real<lower = 0> crcr_covar_trial_corr_eta;