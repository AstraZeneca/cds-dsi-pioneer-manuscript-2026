// Relocated main model (was tumor/sf-ssm-log-space.stan)
functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "hierarchy.stanfunctions"
  #include "full_model.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "_burden.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
  #include "modules/laplace_surrogate/surrogate.stanfunctions"
}

data {
  #include "_hierarchy_data.stan"
  #include "_visit_data.stan"
  #include "_full_model_data.stan"
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
  #include "modules/laplace_surrogate/flags.stan"
  #include "modules/laplace_surrogate/data.stan"

  int<lower = 0, upper = 1> fit_multistate_data;
}

transformed data {
  #include "_hierarchy_transformed_data.stan"
  #include "_forecast_routing_transformed_data.stan"
  int max_all_t = max(max(t_patient_visits) + 1, extend_max_all_t);
  int<lower=0> max_t_width = max_all_t - min(t_patient_visits) + 1;
  #include "_visit_transformed_data.stan"
  // Tumor model has no inline TV-covariate populator (only the dense
  // states_full_grid path is implemented in _ms_burden_tv_covar.stan).
  int has_inline_tv_covar = 0;
  #include "_full_model_transformed_data.stan"
  #include "modules/visits/transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "modules/laplace_surrogate/transformed_data.stan"
  #include "_tumor_observed_covar_transformed_data.stan"
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
  // Burden interface: tumor SLD as the generic burden marker.
  // (median_log_burden_obs / iqr_log_burden_obs are declared in transformed
  // data via _tumor_observed_covar_transformed_data.stan — already in scope.)
  vector[n_patients] log_baseline_burden = log_baseline_sld;
  #include "_ms_burden_tv_covar.stan"
  #include "modules/multistate/transformed_parameters.stan"
  #include "_ms_burden_inline_tv_covar.stan"
  #include "modules/multistate/cond_surv_transform.stan"
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
        for (j in 1:n_forecast_patients) {
          int p = forecast_patient_idx[j];
          int data_start, data_end;
          (data_start, data_end) = get_pos(patient_visit_pos, p);
          int state_start, state_end;
          (state_start, state_end) = get_pos(forecast_visit_pos, j);
          normalized_sld[data_start:data_end] ~ sf_log_space_obs(states[state_start:state_end], measure_sd_sld, log_lod - log_baseline_sld[p]);
        }
      }
    }

    if (fit_multistate_data) {
      profile("multistate loglik") {
        // Slice all per-patient arrays to forecast patients only:
        // log_cond_surv_* are already [n_forecast_patients, max_t];
        // ms_time_*, ms_censored_*, ms_final_state, etc. are [n_patients] and must be sliced.
        ms_final_state[forecast_patient_idx] ~ multistate(
          ones_vector(n_forecast_patients),
          enable_ms_01, enable_ms_02, enable_ms_12, ms_time_scale_12,
          enable_ms_03, enable_ms_32,
          ms_time_01[forecast_patient_idx], ms_time_02[forecast_patient_idx], ms_time_12[forecast_patient_idx],
          ms_time_03[forecast_patient_idx], ms_time_32[forecast_patient_idx],
          ms_censored_01[forecast_patient_idx],
          ms_prog_deterministic[forecast_patient_idx],
          ms_ic_gap_01[forecast_patient_idx],
          t_patient_visits,
          patient_visit_pos,
          log_cond_surv_01,
          log_cond_surv_02,
          log_cond_surv_12_s,
          log_cond_surv_12_t,
          log_cond_surv_03,
          log_cond_surv_32,
          enable_ms_visit_gated_01
        );
      }
    }

    #include "modules/laplace_surrogate/likelihood.stan"
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
  vector[n_forecast_total_groups] tr_effect_level_intercept;
  for (lv in 1:n_levels) {
    int lv_start_output = n_forecast_level_pos[lv];
    int lv_end_output = n_forecast_level_pos[lv + 1] - 1;
    if (enable_level_intercept_tr[lv]) {
      // Index into compacted parameter array using enabled position array
      int lv_start_param = enabled_level_pos_tr_intercept[lv];
      int lv_end_param = enabled_level_pos_tr_intercept[lv + 1] - 1;
      tr_effect_level_intercept[lv_start_output:lv_end_output] =
        tr_sd_level_intercept[lv] * tr_raw_level_intercept[lv_start_param:lv_end_param];
    } else {
      tr_effect_level_intercept[lv_start_output:lv_end_output] =
        rep_vector(0, n_forecast_groups_per_level[lv]);
    }
  }

  // Log rates for all groups at all levels (flattened structure)
  // Each group's rate = population rate + that group's intercept effect
  vector[n_forecast_total_groups] level_log_total_rate = tr_loc_pop + tr_effect_level_intercept;
  vector[n_forecast_total_groups] level_log_decrease_rate = level_log_total_rate + pop_log_decrease_frac;
  vector[n_forecast_total_groups] level_log_growth_rate = level_log_total_rate + pop_log_growth_frac;

  // Residuals for all groups at all levels (vs population)
  vector[n_forecast_total_groups] level_log_growth_rate_residual = level_log_growth_rate - pop_log_growth_rate;
  vector[n_forecast_total_groups] level_log_decrease_rate_residual = level_log_decrease_rate - pop_log_decrease_rate;

  // Patient-level residuals: forecast patients only (background patients have no explicit states)
  matrix[n_forecast_patients, max_t_width] patient_log_growth_rate_residual;
  matrix[n_forecast_patients, max_t_width] patient_log_decrease_rate_residual;

  {
    // Get parent level rates for each forecast patient
    vector[n_forecast_patients] parent_log_growth_rate;
    vector[n_forecast_patients] parent_log_decrease_rate;

    if (n_levels > 1) {
      // Parent is level n_levels - 1
      int parent_lv = n_levels - 1;
      int parent_lv_start, parent_lv_end;
      (parent_lv_start, parent_lv_end) = get_pos(n_forecast_level_pos, parent_lv);

      // Extract parent level rates, then index by forecast patient's group membership
      vector[n_forecast_groups_per_level[parent_lv]] parent_level_growth = level_log_growth_rate[parent_lv_start:parent_lv_end];
      vector[n_forecast_groups_per_level[parent_lv]] parent_level_decrease = level_log_decrease_rate[parent_lv_start:parent_lv_end];
      parent_log_growth_rate = parent_level_growth[patient_level_groups[forecast_patient_idx, parent_lv]];
      parent_log_decrease_rate = parent_level_decrease[patient_level_groups[forecast_patient_idx, parent_lv]];
    } else {
      // No intermediate levels, compare to population
      parent_log_growth_rate = rep_vector(pop_log_growth_rate, n_forecast_patients);
      parent_log_decrease_rate = rep_vector(pop_log_decrease_rate, n_forecast_patients);
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

  // Set generic measure_sd for state_space module
  real measure_sd_obs = measure_sd_sld;

  // Biomarker-agnostic trajectory generation
  #include "modules/state_space/generated_quantities.stan"

  // SLD-specific: aliases + RECIST categorization at assessment visits
  #include "_tumor_categorization.stan"

  // Endpoint output declarations (GQ scope — saved to CSV)
  array[n_forecast_patients] int<lower = 0> sample_target_pfs, spop_target_pfs, sample_ms_pfs, spop_ms_pfs,
                                       spop_target_obs_cens_pfs, sample_pfs, spop_pfs;
  array[n_forecast_patients] int<lower = 0, upper = 1>
    sample_target_right_censored, spop_target_right_censored, spop_target_obs_cens_right_censored,
    sample_ms_right_censored, spop_ms_right_censored,
    sample_right_censored, spop_right_censored;
  array[n_forecast_patients] int<lower = 0> sample_os, spop_os;
  array[n_forecast_patients] int<lower = 0, upper = 1> sample_os_censored, spop_os_censored;
  array[n_forecast_patients] int<lower = 0, upper = 1> spop_is_dropout, sample_is_dropout;
  array[n_forecast_patients] int<lower = 0> spop_dropout_week;
  array[sum(target_right_censored[forecast_patient_idx])] int<lower = 0> forecast_target_pfs;
  array[sum(target_right_censored[forecast_patient_idx])] int<lower = 0, upper = 1> forecast_target_right_censored;
  array[n_forecast_patients] int<lower = 0, upper = 1> sample_target_confirmed_response, spop_target_confirmed_response;
  array[n_forecast_patients] int<lower = 0, upper = 1> sample_target_unconfirmed_response, spop_target_unconfirmed_response;
  vector<lower = 0, upper = 1>[n_trials] sample_target_orr, spop_target_orr;
  vector<lower = 0, upper = 1>[n_cond_group] cond_sample_target_orr = zeros_vector(n_cond_group),
                                              cond_spop_target_orr = zeros_vector(n_cond_group);
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_target_km_est, spop_target_km_est,
                                                               spop_target_obs_cens_km_est,
                                                               sample_ms_pfs_km_est, spop_ms_pfs_km_est,
                                                               sample_pfs_km_est, spop_pfs_km_est;
  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_target_km_est, cond_spop_target_km_est,
                                                                   cond_spop_target_obs_cens_km_est,
                                                                   cond_sample_ms_pfs_km_est, cond_spop_ms_pfs_km_est,
                                                                   cond_sample_pfs_km_est, cond_spop_pfs_km_est;
  array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_target_pfs_n, spop_target_pfs_n,
                                                                  sample_ms_pfs_n, spop_ms_pfs_n,
                                                                  sample_pfs_n, spop_pfs_n;
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_target_pfs_n, cond_spop_target_pfs_n,
                                                                      cond_sample_ms_pfs_n, cond_spop_ms_pfs_n,
                                                                      cond_sample_pfs_n, cond_spop_pfs_n;
  array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_target_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      spop_target_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      sample_ms_pfs_quant     = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      spop_ms_pfs_quant       = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      sample_pfs_quant        = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      spop_pfs_quant          = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
  array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_target_pfs_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_spop_target_pfs_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_sample_ms_pfs_quant     = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_spop_ms_pfs_quant       = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_sample_pfs_quant        = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_spop_pfs_quant          = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
  array[n_trials, n_pfs_quantiles] int sample_target_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       spop_target_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       sample_ms_pfs_quant_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       spop_ms_pfs_quant_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       sample_pfs_quant_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       spop_pfs_quant_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
  array[n_cond_group, n_pfs_quantiles] int cond_sample_target_pfs_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_spop_target_pfs_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_sample_ms_pfs_quant_exceeds_max     = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_spop_ms_pfs_quant_exceeds_max       = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_sample_pfs_quant_exceeds_max        = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_spop_pfs_quant_exceeds_max          = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1] sample_os_km_est, spop_os_km_est;
  array[n_cond_group] vector<lower = 0, upper = 1>[max_all_t + 1] cond_sample_os_km_est, cond_spop_os_km_est;
  array[n_trials] vector<lower = 0>[n_pfs_quantiles] sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_trials),
                                                      spop_os_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_trials);
  array[n_trials, n_pfs_quantiles] int sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_trials),
                                       spop_os_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_trials);
  array[n_cond_group] vector<lower = 0>[n_pfs_quantiles] cond_sample_os_quant = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group),
                                                          cond_spop_os_quant   = rep_array(zeros_vector(n_pfs_quantiles), n_cond_group);
  array[n_cond_group, n_pfs_quantiles] int cond_sample_os_quant_exceeds_max = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group),
                                           cond_spop_os_quant_exceeds_max   = rep_array(zeros_int_array(n_pfs_quantiles), n_cond_group);
  array[n_trials] vector<lower = 0, upper = 1>[n_pfs_timepoints] sample_os_n, spop_os_n;
  array[n_cond_group] vector<lower = 0, upper = 1>[n_pfs_timepoints] cond_sample_os_n, cond_spop_os_n;
  array[n_trials] vector<lower = 0, upper = 1>[max_all_t + 1]
    spop_cif_01   = rep_array(zeros_vector(max_all_t + 1), n_trials),
    spop_cif_02   = rep_array(zeros_vector(max_all_t + 1), n_trials),
    spop_cif_03   = rep_array(zeros_vector(max_all_t + 1), n_trials),
    sample_cif_01 = rep_array(zeros_vector(max_all_t + 1), n_trials),
    sample_cif_02 = rep_array(zeros_vector(max_all_t + 1), n_trials),
    sample_cif_03 = rep_array(zeros_vector(max_all_t + 1), n_trials);

  // Shared endpoint computation (contract vars are local — not saved to CSV)
  {
    array[n_total_visits] int obs_biomarker_cat = recist;
    array[n_total_visits] int rep_biomarker_cat = rep_recist;
    array[n_total_forecast_obs_visits] int forecast_obs_biomarker_cat = forecast_obs_recist;
    array[n_patients] int burden_pfs = pfs;
    array[n_patients] int burden_right_censored = right_censored;
    array[n_patients] int burden_target_pfs = target_pfs;
    array[n_patients] int burden_target_right_censored = target_right_censored;
    // Tumor model does not support visit-gated 0->1 (no PSA covariate)
    int burden_enable_ms_visit_gated_01 = 0;
    real burden_tv_coef_01_val = 0.0;
    vector[0] burden_forecast_obs_log_psa;
    real burden_median_log_psa_obs = 0.0;
    real burden_iqr_log_psa_obs = 0.0;
    // RECIST ORR = confirmed response (>=2 assessments at PR/CR)
    int burden_orr_use_confirmed_response = 1;

    #include "modules/state_space/burden_endpoints.stan"
    #include "modules/state_space/_trial_aggregate_metrics.stan"
  }

  // RECIST accuracy metrics
  #include "modules/tumor/generated_quantities.stan"

  // Back-transformed covariate coefficients (population + arm-level)
  #include "modules/multistate/covar_coef_generated_quantities.stan"
}
