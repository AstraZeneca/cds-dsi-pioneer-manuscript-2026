
row_vector[n_causes] log_lambda_gp_pop_intercept;
vector<lower = 0>[n_causes] log_lambda_gp_pop_alpha;
vector<lower = 0>[n_causes] log_lambda_gp_pop_rho;
array[n_causes] row_vector[max_all_t] log_lambda_gp_pop_eta;

// Trial level hierarchical effect on baseline hazard 
vector<lower = 0>[add_trial_level_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_alpha;
vector<lower = 0>[add_trial_level_baseline_hazard ? n_trials : 0] log_lambda_gp_trial_rho;
array[n_causes] matrix[add_trial_level_baseline_hazard ? n_trials : 0, max_all_t] log_lambda_gp_trial_eta;
array[n_causes] vector[add_trial_level_baseline_hazard ? n_trials : 0] raw_log_lambda_gp_trial_intercept;
vector<lower = 0>[add_trial_level_baseline_hazard ? n_causes : 0] log_lambda_gp_trial_intercept_sd;

