// Relocated main model (was tumor/sf-ssm-log-space.stan)
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "modules/state_space/functions.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "_base_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/state_space/data.stan"
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
  #include "_base_transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/measurement/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"
  #include "modules/other_events/transformed_data.stan"
  #include "modules/state_space/checks.stan"
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
  #include "modules/state_space/transformed_parameters.stan"
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
      profile("tumor loglik") {
        for (i in 1:n_patients) {
          int visit_start, visit_end;
          (visit_start, visit_end) = get_pos(patient_visit_pos, i);
          normalized_sld[visit_start:visit_end] ~ sf_log_space_obs(states[visit_start:visit_end], measure_sd, log_lod - log_baseline_sld[i]);
        }
      }
    }

    if (fit_other_events_data) {
      profile("other events loglik") {
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
}

generated quantities {
  real pop_log_decrease_frac = log_inv_logit(frac_logit_loc_pop);
  real pop_log_growth_frac   = log1m_inv_logit(frac_logit_loc_pop);
  real pop_log_decrease_rate = tr_loc_pop + pop_log_decrease_frac;
  real pop_log_growth_rate   = tr_loc_pop + pop_log_growth_frac;

  // Compute scaled intercept effects for all levels (flattened structure)
  // Note: tr_raw_level_intercept is sized by enabled groups only, so we use
  // enabled_level_pos_tr_intercept for indexing into it
  vector[n_total_groups] tr_effect_level_intercept;
  for (lv in 1:n_levels) {
    int lv_start_output = level_pos[lv];
    int lv_end_output = level_pos[lv + 1] - 1;
    if (enable_level_intercept_tr[lv]) {
      // Index into compacted parameter array using enabled position array
      int lv_start_param = enabled_level_pos_tr_intercept[lv];
      int lv_end_param = enabled_level_pos_tr_intercept[lv + 1] - 1;
      tr_effect_level_intercept[lv_start_output:lv_end_output] =
        tr_sd_level_intercept[lv] * tr_raw_level_intercept[lv_start_param:lv_end_param];
    } else {
      tr_effect_level_intercept[lv_start_output:lv_end_output] =
        rep_vector(0, n_groups_per_level[lv]);
    }
  }

  // Log rates for all groups at all levels (flattened structure)
  // Each group's rate = population rate + that group's intercept effect
  vector[n_total_groups] level_log_total_rate = tr_loc_pop + tr_effect_level_intercept;
  vector[n_total_groups] level_log_decrease_rate = level_log_total_rate + pop_log_decrease_frac;
  vector[n_total_groups] level_log_growth_rate = level_log_total_rate + pop_log_growth_frac;

  // Residuals for all groups at all levels (vs population)
  vector[n_total_groups] level_log_growth_rate_residual = level_log_growth_rate - pop_log_growth_rate;
  vector[n_total_groups] level_log_decrease_rate_residual = level_log_decrease_rate - pop_log_decrease_rate;

  // Patient-level residuals: compare to parent level (level n_levels - 1, or population if n_levels == 1)
  matrix[n_patients, max_t_width] patient_log_growth_rate_residual;
  matrix[n_patients, max_t_width] patient_log_decrease_rate_residual;

  {
    // Get parent level rates for each patient
    vector[n_patients] parent_log_growth_rate;
    vector[n_patients] parent_log_decrease_rate;

    if (n_levels > 1) {
      // Parent is level n_levels - 1
      int parent_lv = n_levels - 1;
      int parent_lv_start, parent_lv_end;
      (parent_lv_start, parent_lv_end) = get_pos(level_pos, parent_lv);

      // Extract parent level rates, then index by patient's group membership
      vector[n_groups_per_level[parent_lv]] parent_level_growth = level_log_growth_rate[parent_lv_start:parent_lv_end];
      vector[n_groups_per_level[parent_lv]] parent_level_decrease = level_log_decrease_rate[parent_lv_start:parent_lv_end];
      parent_log_growth_rate = parent_level_growth[patient_level_groups[, parent_lv]];
      parent_log_decrease_rate = parent_level_decrease[patient_level_groups[, parent_lv]];
    } else {
      // No intermediate levels, compare to population
      parent_log_growth_rate = rep_vector(pop_log_growth_rate, n_patients);
      parent_log_decrease_rate = rep_vector(pop_log_decrease_rate, n_patients);
    }

    if (enable_patient_process_noise_tr) {
      // Time-varying rates: direct subtraction
      patient_log_growth_rate_residual = patient_log_growth_rate - rep_matrix(parent_log_growth_rate, max_t_width);
      patient_log_decrease_rate_residual = patient_log_decrease_rate - rep_matrix(parent_log_decrease_rate, max_t_width);
    } else {
      // Constant rates: broadcast single column across all time points
      patient_log_growth_rate_residual = patient_log_growth_rate[, 1] * ones_row_vector(max_t_width) - rep_matrix(parent_log_growth_rate, max_t_width);
      patient_log_decrease_rate_residual = patient_log_decrease_rate[, 1] * ones_row_vector(max_t_width) - rep_matrix(parent_log_decrease_rate, max_t_width);
    }
  }

  #include "_endpoints_generated_quantities.stan"  
  #include "modules/state_space/generated_quantities.stan"
}
