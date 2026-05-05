// Relocated main model (was tumor/sf-ssm-log-space.stan)
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "_base_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/visits/data.stan"
  #include "modules/tumor/hyperparams.stan"
  #include "modules/state_space/data.stan"
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  int<lower = 0, upper = 1> fit_multistate_data;
}

transformed data {
  #include "_base_transformed_data.stan"
  #include "modules/visits/transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "modules/state_space/checks.stan"
}

parameters {
  #include "modules/tumor/parameters.stan"
  #include "modules/multistate/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "modules/state_space/transformed_parameters.stan"
  #include "_ms_time_varying_covar.stan"
  #include "modules/multistate/transformed_parameters.stan"
}

model {
  #include "modules/tumor/priors.stan"
  #include "modules/multistate/priors.stan"
  #include "modules/tr/priors.stan"
  #include "modules/frac/priors.stan"
  #include "modules/init/priors.stan"

  profile("loglik") {
    if (fit_tumor_data) {
      profile("tumor loglik") {
        for (i in 1:n_patients) {
          int visit_start, visit_end;
          (visit_start, visit_end) = get_pos(patient_visit_pos, i);
          normalized_sld[visit_start:visit_end] ~ sf_log_space_obs(states[visit_start:visit_end], measure_sd_sld, log_lod - log_baseline_sld[i], measure_nu_sld);
        }
      }
    }

    if (fit_multistate_data) {
      profile("multistate loglik") {
        ms_final_state ~ multistate(
          enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
          enable_ms_03, enable_ms_32,
          ms_time_01, ms_time_02, ms_time_12,
          ms_time_03, ms_time_32,
          ms_censored_01, ms_censored_02, ms_censored_12,
          ms_censored_32,
          ms_prog_deterministic,
          ms_ic_gap_01,
          t_patient_visits,
          patient_visit_pos,
          log_cond_surv_01,
          log_cond_surv_02,
          log_cond_surv_12_s,
          log_cond_surv_12_t,
          log_cond_surv_03,
          log_cond_surv_32
        );
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
