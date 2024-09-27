// Baseline hazard GP parameters
vector<lower = 0>[n_causes] log_crcr_lambda_gp_alpha;
vector<lower = 0>[n_causes] log_crcr_lambda_gp_rho;
matrix[max_confresp_week, n_causes] log_crcr_lambda_gp_eta;
vector[n_causes] log_crcr_lambda_gp_intercept;

// Trial level hierarchical effect on baseline hazard 
vector<lower = 0>[add_trial_level ? n_causes : 0] log_crcr_lambda_gp_trial_alpha;
vector<lower = 0>[add_trial_level ? n_causes : 0] log_crcr_lambda_gp_trial_rho;
array[add_trial_level ? n_trials : 0] matrix[max_confresp_week, n_causes] log_crcr_lambda_gp_trial_eta;
matrix[add_trial_level ? n_trials : 0, n_causes] raw_log_crcr_lambda_gp_trial_intercept;
vector[add_trial_level ? n_causes : 0] log_crcr_lambda_gp_trial_intercept_sd;

// GLM parameters
matrix[n_tumor_covar, n_causes] crcr_tumor_stim_pop_coef;
matrix[n_covar, n_causes] crcr_covar_effect;  

// GLM hierarchical parameters
vector<lower = 0>[add_trial_level && add_trial_level_glm ? n_tumor_covar + n_covar : 0] crcr_covar_trial_sd;
cholesky_factor_corr[add_trial_level ? n_tumor_covar + n_covar : 0] L_crcr_covar_trial_corr;
array[add_trial_level && add_trial_level_glm ? n_trials : 0] matrix[n_tumor_covar + n_covar, n_causes] raw_crcr_covar_trial_coef;