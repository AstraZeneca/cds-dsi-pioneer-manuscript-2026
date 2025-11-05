// Relocated main model (was tumor/sf-ssm-log-space.stan)
functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "_sf_functions.stan"
  #include "../recist.stanfunctions"
}

data {
  #include "../base_data.stan"
  #include "../tumor/base_data.stan" // tumor-specific base
  #include "_sf_outcomes_info.stan"
  #include "legacy/sf-ssls-hyperparam.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"
}

transformed data {
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "_sf_transformed_data.stan"
  #include "_other_events_transformed_data.stan"
  #include "_sf-checks.stan"
}

parameters {
  #include "_other_events_parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
  #include "legacy/sf-ssls-parameters.stan"
}

transformed parameters {
  #include "_other_events_transformed_parameters.stan"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "legacy/sf-ssls-transformed_parameters.stan"

  matrix[n_patients, n_causes] patient_response_lp = rep_matrix(0, n_patients, n_causes); 
  patient_response_lp[, 1] = calc_pch_loglik(
    other_events_pfs, 
    other_events_right_censored, 
    other_events_interval_censored,
    0, 
    log_cond_prob_surv[1]
  );
}

model {
  #include "_other_events_priors.stan"
  #include "modules/tr/priors.stan"
  #include "modules/frac/priors.stan"
  #include "modules/init/priors.stan"
  #include "legacy/sf-ssls-priors.stan"

  profile("loglik") { 
    if (fit_tumor_data) {
      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        normalized_sld[visit_start:visit_end] ~ sf_log_space_obs(states[visit_start:visit_end], measure_sd, log_lod - log(sum_tumor_size[visit_start]));
      }

      for (s in 1:n_trials) for (k in 1:n_causes) target += sum(get_sub_vector(patient_response_lp[, k], trial_patient_pos, s));
    }
  }
}

generated quantities {
  real pop_log_decrease_frac = log_inv_logit(frac_logit_loc_pop);
  real pop_log_growth_frac   = log1m_inv_logit(frac_logit_loc_pop);
  real pop_log_decrease_rate = tr_loc_pop + pop_log_decrease_frac;
  real pop_log_growth_rate   = tr_loc_pop + pop_log_growth_frac;
  vector[n_trials] trial_log_total_rate = tr_loc_pop + tr_effect_trial_intercept;
  vector[n_trials] trial_log_decrease_rate = trial_log_total_rate + pop_log_decrease_frac;
  vector[n_trials] trial_log_growth_rate   = trial_log_total_rate + pop_log_growth_frac;
  vector[n_trials] trial_log_growth_rate_residual = trial_log_growth_rate - pop_log_growth_rate;
  vector[n_patients] patient_log_growth_rate_residual = patient_log_growth_rate - trial_log_growth_rate[patient_trial];
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  vector[n_patients] patient_log_decrease_rate_residual = patient_log_decrease_rate - trial_log_decrease_rate[patient_trial];

  vector<lower = 0, upper = 1>[max_all_t] all_growth_factor = get_growth_lag_factor(all_tumor_measure_t, exp(pop_log_growth_lag), exp(pop_log_growth_transition_rate));
  
  matrix[max_all_t, 2] all_scaled_process_sd = scale_process_sd(all_tumor_measure_t, pop_process_sd);

  #include "_endpoints_generated_quantities.stan"  
  #include "_sf_accuracy_generated_quantities.stan"
}
