vector[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_intercept_mean;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_intercept_sd;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_alpha_sd;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_rho_alpha;
vector<lower = 0>[separate_baseline_hazard ? n_trials : 1] log_lambda_gp_rho_beta;

real<lower = 0> log_lambda_gp_trial_alpha_sd;
real<lower = 0> log_lambda_gp_trial_intercept_sd_sd;
  