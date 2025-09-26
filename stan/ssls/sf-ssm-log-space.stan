// Relocated main model (was tumor/sf-ssm-log-space.stan)
functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "legacy/sf-ssls_functions.stan"
  #include "../recist.stanfunctions"
}

data {
  #include "../base_data.stan"
  #include "../tumor/base_data.stan" // tumor-specific base
  #include "legacy/sf-ssls-outcomes_info.stan"
  #include "legacy/sf-ssls-hyperparam.stan"
  // Module hyperparameter includes (new convention). For now these
  // names are still satisfied via aliasing in transformed data until
  // R side switches to providing them directly.
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  // Total rate flags
  int<lower=0,upper=1> enable_pop_cov_tr;
  int<lower=0,upper=1> enable_trial_intercept_tr;
  int<lower=0,upper=1> enable_trial_cov_tr;
  int<lower=0,upper=1> enable_patient_intercept_tr;
  int<lower=0,upper=1> enable_patient_cov_tr;
  // Fraction flags
  int<lower=0,upper=1> enable_pop_cov_frac;
  int<lower=0,upper=1> enable_trial_intercept_frac;
  int<lower=0,upper=1> enable_trial_cov_frac;
  int<lower=0,upper=1> enable_patient_intercept_frac;
  int<lower=0,upper=1> enable_patient_cov_frac;
  // Initial proportion flags
  int<lower=0,upper=1> enable_pop_cov_init;
  int<lower=0,upper=1> enable_trial_intercept_init;
  int<lower=0,upper=1> enable_trial_cov_init;
  int<lower=0,upper=1> enable_patient_intercept_init;
  int<lower=0,upper=1> enable_patient_cov_init;
}

transformed data {
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "_sf_transformed_data.inc"
  #include "_other_events_transformed_data.inc"
  #include "_mature_cutoffs_transformed_data.inc"
  #include "legacy/sf-ssls-outcomes_info_transformed_data.stan"
  #include "sf-checks.stan"
}

parameters {
  #include "_other_events_parameters.inc"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
  #include "legacy/sf-ssls-parameters.stan"
}

transformed parameters {
  #include "_other_events_transformed_parameters.inc"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "legacy/sf-ssls-transformed_parameters.stan"

  matrix[n_train_patients, n_causes] patient_response_lp = rep_matrix(0, n_train_patients, n_causes); 
  patient_response_lp[, 1] = calc_pch_loglik(
    non_target_pfs[train_patients_pos:train_patients_end], 
    non_target_right_censored[train_patients_pos:train_patients_end], 
    zeros_int_array(n_train_patients),
    0, 
    log_cond_prob_surv[1]
  );
}

model {
  #include "_other_events_priors.inc"
  #include "modules/tr/priors.stan"
  #include "modules/frac/priors.stan"
  #include "modules/init/priors.stan"
  #include "legacy/sf-ssls-priors.stan"

  profile("loglik") { 
    if (fit_tumor_data) {
      for (i in train_patients_pos:train_patients_end) {
        int train_idx = i - train_patients_pos + 1;
        int train_visit_start, train_visit_end;
        (train_visit_start, train_visit_end) = get_pos(train_patient_visit_pos, train_idx);
        int visit_pos, visit_end;
        (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
        normalized_sld[visit_pos:visit_end] ~ sf_log_space_obs(states[train_visit_start:train_visit_end], measure_sd, log_lod - log(sum_tumor_size[visit_pos]));
      }
      for (s in 1:n_trials) for (k in 1:n_causes) target += sum(get_sub_vector(patient_response_lp[, k], train_trial_patient_pos, s));
    }
  }
}

generated quantities {
  real pop_log_decrease_frac = -log1p_exp(-pop_decrease_frac_logit);
  real pop_log_growth_frac   = -log1p_exp(pop_decrease_frac_logit);
  real pop_log_decrease_rate = pop_log_total_rate + pop_log_decrease_frac;
  real pop_log_growth_rate   = pop_log_total_rate + pop_log_growth_frac;
  vector[n_trials] trial_log_total_rate = pop_log_total_rate + trial_log_total_rate_effect;
  vector[n_trials] trial_log_decrease_rate = trial_log_total_rate + pop_log_decrease_frac;
  vector[n_trials] trial_log_growth_rate   = trial_log_total_rate + pop_log_growth_frac;
  vector[n_trials] trial_log_growth_rate_residual = trial_log_growth_rate - pop_log_growth_rate;
  vector[n_train_patients] patient_log_growth_rate_residual = patient_log_growth_rate - trial_log_growth_rate[patient_trial[train_patients_pos:train_patients_end]];
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  vector[n_train_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - trial_log_decrease_rate[patient_trial[train_patients_pos:train_patients_end]];

  // #include "legacy/sf-ssls-accuracy_gen_quant.stan"
}
