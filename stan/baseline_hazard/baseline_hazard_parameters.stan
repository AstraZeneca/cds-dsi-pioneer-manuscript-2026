// Baseline hazard GP parameters
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_alpha;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_rho;
array[separate_baseline_hazard ? n_trials : 1] vector[max_all_t] log_lambda_gp_eta;
vector[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_intercept;

// Trial level hierarchical effect on baseline hazard 
real<lower = 0> log_lambda_gp_trial_alpha;
real<lower = 0> log_lambda_gp_trial_rho;
array[add_trial_level_baseline_hazard ? n_trials : 0] vector[max_all_t] log_lambda_gp_trial_eta;
vector[add_trial_level_baseline_hazard ? n_trials : 0] raw_log_lambda_gp_trial_intercept;
real<lower = 0> log_lambda_gp_trial_intercept_sd;
 