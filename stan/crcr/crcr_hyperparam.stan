real log_crcr_lambda_gp_intercept_mean;
real<lower = 0> log_crcr_lambda_gp_intercept_sd;
real<lower = 0> log_crcr_lambda_gp_alpha_sd;
real<lower = 0> log_crcr_lambda_gp_trial_alpha_sd;
real<lower = 0> log_crcr_lambda_gp_rho_alpha;
real<lower = 0> log_crcr_lambda_gp_rho_beta;
real<lower = 0> log_crcr_lambda_gp_trial_intercept_sd_sd;

vector<lower = 0>[2] crcr_tumor_stim_pop_coef_sd;
vector<lower = 0>[n_covar] crcr_covar_effect_sd;