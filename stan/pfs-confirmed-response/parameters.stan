#include "../baseline_hazard/baseline_hazard_parameters.stan"
#include "../crcr/crcr_parameters.stan"

// Proportional hazard parameters

array[no_prop_hazard ? 0 : n_prop_separate_trials] vector[n_tumor_covar] tumor_stim_pop_coef;
array[no_prop_hazard ? 0 : n_prop_separate_trials] vector[n_covar] covar_effect;  
vector[no_prop_hazard ? 0 : n_prop_separate_trials] conf_resp_effect;

vector<lower = 0>[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar + 1 : 0] covar_trial_sd;
cholesky_factor_corr[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar + 1 : 0] L_covar_trial_corr;
array[add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar + 1] raw_covar_trial_coef;