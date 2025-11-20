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
  #include "modules/measurement/hyperparams.stan"
  #include "modules/other_events/data.stan"
  #include "modules/other_events/hyperparams.stan"
  #include "modules/other_events/flags.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/measurement/flags.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  int<lower = 0, upper = 1> fit_other_events_data;
}

transformed data {
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "modules/measurement/transformed_data.stan"
  #include "_sf_transformed_data.stan"
  #include "modules/other_events/transformed_data.stan"
  #include "_sf-checks.stan"
}

parameters {
  #include "modules/measurement/parameters.stan"
  #include "modules/other_events/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/measurement/transformed_parameters.stan"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "_sf_transformed_parameters.stan"
  #include "modules/other_events/transformed_parameters.stan"
}

model {
  #include "modules/measurement/priors.stan"
  #include "modules/other_events/priors.stan"
  #include "modules/tr/priors.stan"
  #include "modules/frac/priors.stan"
  #include "modules/init/priors.stan"

  profile("loglik") { 
    if (fit_tumor_data) {
      for (i in 1:n_patients) {
        int visit_start, visit_end;
        (visit_start, visit_end) = get_pos(patient_visit_pos, i);
        normalized_sld[visit_start:visit_end] ~ sf_log_space_obs(states[visit_start:visit_end], measure_sd, log_lod - log_baseline_sld[i]);
      }
    }

    if (fit_other_events_data) {
      // Other events likelihood contribution
      matrix[n_patients, n_causes] patient_response_lp = rep_matrix(0, n_patients, n_causes); 
      patient_response_lp[, 1] = calc_pch_loglik(
        ic_other_events_pfs, 
        other_events_right_censored, 
        zeros_int_array(n_patients), // other_events_interval_censored
        0, 
        log_cond_prob_surv[1]
      );

      target += sum(patient_response_lp);
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
  matrix[n_patients, max_t_width] patient_log_growth_rate_residual = 
    patient_log_growth_rate - rep_matrix(trial_log_growth_rate[patient_trial], max_t_width);
  vector[n_trials] trial_log_decrease_rate_residual = trial_log_decrease_rate - pop_log_decrease_rate;
  matrix[n_patients, max_t_width] patient_log_decrease_rate_residual = 
    patient_log_decrease_rate - rep_matrix(trial_log_decrease_rate[patient_trial], max_t_width);

  #include "_endpoints_generated_quantities.stan"  
  #include "_sf_accuracy_generated_quantities.stan"
}
