// Baseline hazard GP parameters
// array[n_base_separate_trials] vector[n_causes] log_crcr_lambda_gp_intercept;
// array[n_base_separate_trials] vector<lower = 0>[n_causes] log_crcr_lambda_gp_alpha;
// array[n_base_separate_trials] vector<lower = 0>[n_causes] log_crcr_lambda_gp_rho;
// array[n_base_separate_trials] matrix[max_confresp_week, n_causes] log_crcr_lambda_gp_eta;
array[n_causes] vector[n_base_separate_trials] log_crcr_lambda_gp_intercept;
array[n_causes] vector<lower = 0>[n_base_separate_trials] log_crcr_lambda_gp_alpha;
array[n_causes] vector<lower = 0>[n_base_separate_trials] log_crcr_lambda_gp_rho;
array[n_causes] matrix[n_base_separate_trials, max_confresp_week] log_crcr_lambda_gp_eta;

// Trial level hierarchical effect on baseline hazard 
vector<lower = 0>[add_trial_level_baseline_hazard ? n_causes : 0] log_crcr_lambda_gp_trial_alpha;
vector<lower = 0>[add_trial_level_baseline_hazard ? n_causes : 0] log_crcr_lambda_gp_trial_rho;
// array[add_trial_level_baseline_hazard ? n_trials : 0] matrix[max_confresp_week, n_causes] log_crcr_lambda_gp_trial_eta;
array[n_causes] matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_confresp_week] log_crcr_lambda_gp_trial_eta;
// matrix[add_trial_level_baseline_hazard ? n_trials : 0, n_causes] raw_log_crcr_lambda_gp_trial_intercept;
array[n_causes] vector[add_trial_level_baseline_hazard ? n_trials : 0] raw_log_crcr_lambda_gp_trial_intercept;
vector<lower = 0>[add_trial_level_baseline_hazard ? n_causes : 0] log_crcr_lambda_gp_trial_intercept_sd;

// GLM parameters
// array[n_prop_separate_trials] matrix[n_tumor_covar, n_causes] crcr_tumor_stim_pop_coef;
// array[n_prop_separate_trials] matrix[n_covar, n_causes] crcr_covar_effect;  
array[n_causes, n_prop_separate_trials] vector[n_tumor_covar] crcr_tumor_stim_pop_coef;
array[n_causes, n_prop_separate_trials] vector[n_covar] crcr_covar_effect;  

// GLM hierarchical parameters
vector<lower = 0>[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar : 0] crcr_covar_trial_sd;
cholesky_factor_corr[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar : 0] L_crcr_covar_trial_corr;
// array[add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] matrix[n_tumor_covar + n_covar, n_causes] raw_crcr_covar_trial_coef;
array[n_causes, add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar] raw_crcr_covar_trial_coef;