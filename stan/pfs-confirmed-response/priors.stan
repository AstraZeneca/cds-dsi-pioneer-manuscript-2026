profile("crcr priors") {
  #include "../crcr/crcr_priors.stan"
}

#include "../baseline_hazard/baseline_hazard_priors.stan"

profile("pfs priors") {
  for (s in 1:n_prop_separate_trials) {
    tumor_stim_pop_coef[s] ~ normal(0, tumor_stim_pop_coef_sd[s]);
    conf_resp_effect[s] ~ normal(conf_resp_effect_mean[s], conf_resp_effect_sd[s]); 
    covar_effect[s] ~ normal(covar_effect_mean[s], covar_effect_sd[s]);
  }
  
  if (add_trial_level_prop_hazard && !no_prop_hazard) {
    covar_trial_sd ~ normal(0, covar_trial_sd_sd);
    L_covar_trial_corr ~ lkj_corr_cholesky(covar_trial_corr_eta);
    
    for (s in 1:n_trials) {
      raw_covar_trial_coef[s] ~ std_normal();
    }
  }
}
