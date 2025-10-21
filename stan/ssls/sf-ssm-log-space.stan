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
  #include "legacy/sf-ssls-outcomes_info_transformed_data.stan"
  #include "sf-checks.stan"
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
    non_target_pfs, 
    non_target_right_censored, 
    zeros_int_array(n_patients),
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
  #include "_endpoints_generated_quantities.stan"  
  #include "legacy/sf-ssls-accuracy_gen_quant.stan"
}
