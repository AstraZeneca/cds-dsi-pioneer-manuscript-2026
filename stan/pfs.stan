functions {
  #include "util.stan"
  #include "pfs_functions.stan"
}

data {
  // Model settings
  int<lower = 0, upper = 1> fit_data; // If 0, just do prior prediction
  int<lower = 0, upper = 1> gen_pfs; // Generate simulated data
  int<lower = 0, upper = 1> gen_interval_censored; // Should the generated PFS be interval censored?
  int<lower = 0, upper = 1> ignore_interval_censoring; // Treat observed PFS as true pfs and ignore t_measure.
  int<lower = 0, upper = 1> use_tumor_model; // Should the tumor model be jointly fit
  int<lower = 0, upper = 1> fit_post_2nd_meaure_only; // Should we exclude all survival intervals before second tumor assessment (post-treatment) from loglik calculation.
  
  // Hierarchical settings 
  int<lower = 0, upper = 1> add_trial_level;
  int<lower = 0, upper = 1> add_tumor_location_level;
  
  // 0: None; 1: First two, no interaction; 2: First two, interaction; 3: First two percentage difference; 4: tumor intercept only 
  // 5: quadratic
  int<lower = 0, upper = 5> tumor_hazard_type; 
 
  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
  
  // This is the tumor model settings for joint modeling. tumor_data is not a good name.
  #include "tumor_data.stan"
    
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;

  // I use these to generate 2d hazard ratios over tumor sizes observed over time. 
  int<lower = 0> n_grid_tumors;
  array[n_grid_tumors] int<lower = 1, upper = sum(n_patient_tumors)> grid_tumors;
  
  // Hyperparam
  #include "baseline_hazard_hyperparam.stan"
  
  real<lower = 0> tumor_stim_pop_intercept_sd;
  vector<lower = 0>[5] tumor_stim_pop_coef_sd;
  vector<lower = 0>[6] tumor_stim_trial_coef_sd_sd;
  vector<lower = 0>[6] tumor_stim_location_coef_sd_sd;
}

transformed data {
  #include "tumor_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  
  if (min(n_patient_screening_t) <= 0) {
    reject("We need at least one screening measurement per patient, for now: ", min(n_patient_screening_t));
  }
    
  int n_covar_col = calc_n_covar_col(tumor_hazard_type); // number of columns in covariates design matrix, excluding the intercept. 
 
  int max_measures = 2; 
  array[n_tumors, max_measures] int<lower = min(t_measure), upper = max(t_measure)> tumor_covar_t; // ts (weeks) of the assessments used in the covar design matrix
  array[n_patients] int<lower = min(t_measure), upper = max(t_measure)> patient_max_2nd_tumor_t; // What t is the second assessment in the covar design matrix 
  matrix[n_tumors, n_covar_col] tumor_covar; // This is the design matrix with the covar in the first two (or _n_) assessments.
  matrix[n_tumors, n_covar_col] uncentered_tumor_covar; // Just scaled
  vector[n_covar_col] tumor_covar_mean;
  vector<lower = 0>[n_covar_col] tumor_covar_sd = rep_vector(0, n_covar_col);
 
  (tumor_covar_t, patient_max_2nd_tumor_t, tumor_covar, uncentered_tumor_covar, tumor_covar_mean, tumor_covar_sd) = 
    prepare_early_tumors_covar(tumor_size, tumor_hazard_type, n_patients, n_tumors, n_patient_tumors, n_measures, t_measure, n_screening_t, max_measures); 
}

parameters {
  // For joint modeling
  #include "tumor_parameters.stan"
  
  #include "baseline_hazard_parameters.stan"
  
  // Per tumor population-level influence on hazard 
  real<lower = 0> tumor_stim_pop_intercept; // DO NOT REMOVE; this is a per tumor intercept and not per patient intercept which is included in lambda.
  row_vector[n_covar_col] tumor_stim_pop_coef;
 
  // Trial level hierarchical tumor effect
  matrix[add_trial_level ? n_trials : 0, n_covar_col + 1] raw_tumor_stim_trial_coef;
  row_vector<lower = 0>[add_trial_level ? n_covar_col + 1 : 0] tumor_stim_trial_coef_sd;
  
  // Organ level hierarchical tumor effect 
  matrix[add_tumor_location_level ? n_tumor_locations : 0, n_covar_col + 1] raw_tumor_stim_location_coef;
  row_vector<lower = 0>[add_tumor_location_level ? n_covar_col + 1 : 0] tumor_stim_location_coef_sd;
}

transformed parameters {
  #include "tumor_transformed_parameters.stan"
  #include "baseline_hazard_transformed_parameters.stan"
  
  // Hazard ratio trial-level parameters 
  vector[n_trials] tumor_stim_trial_intercept = add_trial_level ? raw_tumor_stim_trial_coef[, 1] * tumor_stim_trial_coef_sd[1] : rep_vector(0, n_trials);
  matrix[n_trials, n_covar_col] tumor_stim_trial_coef;
  
  for (s in 1:n_trials) {
    if (n_covar_col > 0) { 
      // TODO Add correlation between intercept and the coefs
      tumor_stim_trial_coef[s] = add_trial_level ? raw_tumor_stim_trial_coef[s, 2:] .* tumor_stim_trial_coef_sd[2:] : rep_row_vector(0, n_covar_col);
    }
    
  }
  
  // Hazard ratio organ-level parameters 
  vector[n_tumor_locations] tumor_stim_location_intercept;
  matrix[n_tumor_locations, n_covar_col] tumor_stim_location_coef;
  
  tumor_stim_location_intercept = add_tumor_location_level ? raw_tumor_stim_location_coef[, 1] * tumor_stim_location_coef_sd[1] : rep_vector(0, n_tumor_locations);
  
  if (n_covar_col > 0) { 
    for (l in 1:n_tumor_locations) {
      // TODO Add correlation between intercept and the coefs
      tumor_stim_location_coef[l] = 
        add_tumor_location_level ? raw_tumor_stim_location_coef[l, 2:] .* tumor_stim_location_coef_sd[2:] : rep_row_vector(0, n_covar_col);
    }
  }
 
  // This is a sum of the contribution of all a patient's tumors to their hazard 
  vector[n_patients] total_time_invar_tumor_stim;
  vector[n_patients] total_time_invar_tumor_stim_no_intercept; // This used outside the model, so don't delete it.
  vector[n_time_periods] disease_progress_pred; // DP predictor for all patients at all observed and unobserved intervals.
  
  (total_time_invar_tumor_stim, total_time_invar_tumor_stim_no_intercept, disease_progress_pred) = calc_disease_progress_pred_from_early_tumors(
    pfs, tumor_covar,
    tumor_hazard_type, n_patients, n_patient_tumors,
    patient_trial,
    tumor_location,
    right_uncensored, interval_censored, 
    log_trial_lambda,
    tumor_stim_pop_intercept, tumor_stim_trial_intercept, tumor_stim_location_intercept, 
    tumor_stim_pop_coef, tumor_stim_trial_coef, tumor_stim_location_coef,
    gen_pfs, max_all_t
  ); 
  
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob = inv_cloglog(disease_progress_pred); // DP conditional probability for all patients (same as above)
  
  // { // Calculate patient-interval conditional probability of disease progression.
  //   int tumor_pos = 1;
  //   int pfs_interval_pos = 1;
  // 
  //   for (i in 1:n_patients) {
  //     int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
  //     int n_intervals = gen_pfs ? max_all_t : pfs[i] + right_uncensored[i] + interval_censored[i];
  //     int pfs_interval_end = pfs_interval_pos + n_intervals - 1; 
  //     array[n_patient_tumors[i]] int patient_tumor_locations = tumor_location[tumor_pos:tumor_end];
  //     
  //     vector[n_patient_tumors[i]] patient_stim_intercept = 
  //       tumor_stim_pop_intercept + tumor_stim_trial_intercept[patient_trial[i]] + tumor_stim_location_intercept[patient_tumor_locations];
  //       
  //     matrix[n_patient_tumors[i], n_covar_col] patient_stim_coef =
  //       // rep_matrix(tumor_stim_pop_coef .* exp(tumor_stim_trial_coef[patient_trial[i]]), n_patient_tumors[i]) .* exp(tumor_stim_location_coef[patient_tumor_locations]);
  //       rep_matrix(tumor_stim_pop_coef + tumor_stim_trial_coef[patient_trial[i]], n_patient_tumors[i]) +tumor_stim_location_coef[patient_tumor_locations];
  //      
  //     vector[n_patient_tumors[i]] tumor_stim = linear_tumor_stimulus(patient_stim_intercept, patient_stim_coef, tumor_covar[tumor_pos:tumor_end]);
  //     total_time_invar_tumor_stim[i] = tumor_hazard_type > 0 ? sum(tumor_stim) : 0;
  //     total_time_invar_tumor_stim_no_intercept[i] = 
  //       tumor_hazard_type > 0 ? sum(linear_tumor_stimulus(rep_vector(0, n_patient_tumors[i]), patient_stim_coef, tumor_covar[tumor_pos:tumor_end])) : 0;
  //    
  //     disease_progress_pred[pfs_interval_pos:pfs_interval_end] = log_trial_lambda[patient_trial[i], 1:n_intervals] + total_time_invar_tumor_stim[i];
  //     
  //     tumor_pos = tumor_end + 1;
  //     pfs_interval_pos = pfs_interval_end + 1;
  //   }
  // }
  // 
  // disease_progress_prob = inv_cloglog(disease_progress_pred); 
}

model {
  // For joint modeling
  #include "tumor_model.stan"
  
  // Priors
  
  #include "baseline_hazard_priors.stan"
  
  tumor_stim_pop_intercept ~ normal(0, tumor_stim_pop_intercept_sd);
 
  if (add_trial_level) { 
    tumor_stim_trial_coef_sd ~ normal(0, tumor_stim_trial_coef_sd_sd[:(n_covar_col + 1)]);
    
    to_vector(raw_tumor_stim_trial_coef) ~ std_normal();
  }
  
  if (add_tumor_location_level) { 
    tumor_stim_location_coef_sd ~ normal(0, tumor_stim_location_coef_sd_sd[:(n_covar_col + 1)]);
    to_vector(raw_tumor_stim_location_coef) ~ std_normal();
  }
  
  if ((tumor_hazard_type > 0 && tumor_hazard_type < 4) || tumor_hazard_type > 4) {
    tumor_stim_pop_coef[1] ~ normal(0, tumor_stim_pop_coef_sd[1]);

    if (tumor_hazard_type == 1 || tumor_hazard_type == 2) {
      tumor_stim_pop_coef[2] ~ normal(0, tumor_stim_pop_coef_sd[2]);

      if (tumor_hazard_type == 2 || tumor_hazard_type == 5) {
        tumor_stim_pop_coef[3] ~ normal(0, tumor_stim_pop_coef_sd[3]);
      }
      
      if (tumor_hazard_type == 5) {
        tumor_stim_pop_coef[4] ~ normal(0, tumor_stim_pop_coef_sd[4]);
        tumor_stim_pop_coef[5] ~ normal(0, tumor_stim_pop_coef_sd[5]);
      }
    }
  }
  
  // Likelihood 
  
  if (fit_data) {
    pfs ~ pch(
      right_uncensored, interval_censored, ignore_interval_censoring, 
      disease_progress_prob, 
      gen_pfs ? max_all_t : 0, 
      fit_post_2nd_meaure_only ? patient_max_2nd_tumor_t : rep_array(1, n_patients) 
    );
  }
}

generated quantities {
  #include "pfs_generated_quant.stan"
  
  vector<lower = 0, upper = 1>[max_all_t] base_pf_cond_prob = 1 - inv_cloglog(log_lambda); // Progression free conditional prob if not using covar
    // Progress free conditional probability if only 1 tumor per patient fixed at size = 1 
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_pf_cond_prob; 
  
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_base_pf_cond_prob; 
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_one_tumor_pf_cond_prob; 
  
  if (add_trial_level) {
    for (s in 1:n_trials) {
      trial_base_pf_cond_prob[s] = 1 - inv_cloglog(log_trial_lambda[s]);
      
      if (tumor_hazard_type > 0) {
        trial_one_tumor_pf_cond_prob[s] = 
          1 - calculate_progress_linear_prob(log_trial_lambda[s], tumor_stim_trial_intercept[s], tumor_stim_trial_coef[s], rep_row_vector(1, n_covar_col));
      } else {
        trial_one_tumor_pf_cond_prob[s] = trial_base_pf_cond_prob[s];
      } 
    }
  }
 
  if (tumor_hazard_type > 0) {
    one_tumor_pf_cond_prob = 
      1 - calculate_progress_linear_prob(log_lambda, tumor_stim_pop_intercept, tumor_stim_pop_coef, rep_row_vector(1, n_covar_col)); 
  } else {
    one_tumor_pf_cond_prob = base_pf_cond_prob;
  } 
  
  vector<lower = 0, upper = 1>[max_all_t] base_survival;
  vector<lower = 0, upper = 1>[max_all_t] one_tumor_survival;
  real<lower = 0, upper = max_all_t> base_cond_expected_pfs;
  real<lower = 0, upper = max_all_t> one_tumor_cond_expected_pfs;
  real<lower = 0, upper = max_all_t> base_cond_median_pfs;
  real<lower = 0, upper = max_all_t> one_tumor_cond_median_pfs;
  
  row_vector<lower = 0, upper = 1>[max_all_t] base_dp_prob;
  row_vector<lower = 0, upper = 1>[max_all_t] one_tumor_dp_prob;
  
  row_vector<lower = 0, upper = 1>[max_all_t] base_dp_prob_not_censored;
  row_vector<lower = 0, upper = 1>[max_all_t] one_tumor_dp_prob_not_censored;
  
  // The below is to evaluate the tumor contribution for the "grid" tumor sizes requested.
  
  array[add_trial_level ? n_trials : 0] vector[n_grid_tumors] trial_tumor_effect;
  array[add_tumor_location_level ? n_tumor_locations : 0] vector[n_grid_tumors] location_tumor_effect;
  vector[n_grid_tumors] pop_tumor_effect = linear_tumor_stimulus(
    // rep_vector(tumor_stim_pop_intercept, n_grid_tumors), 
    rep_vector(0, n_grid_tumors), 
    rep_matrix(tumor_stim_pop_coef, n_grid_tumors), 
    tumor_covar[grid_tumors]
  );
  
  if (add_trial_level) { 
    for (s in 1:n_trials) {
      trial_tumor_effect[s] = linear_tumor_stimulus(
        // rep_vector(tumor_stim_pop_intercept + tumor_stim_trial_intercept[s], n_grid_tumors), 
        rep_vector(0, n_grid_tumors), 
        // rep_matrix(tumor_stim_pop_coef .* exp(tumor_stim_trial_coef[s]), n_grid_tumors), 
        rep_matrix(tumor_stim_pop_coef + tumor_stim_trial_coef[s], n_grid_tumors), 
        tumor_covar[grid_tumors]
      );
    }
  }
    
  if (add_tumor_location_level) { 
    for (l in 1:n_tumor_locations) {
      location_tumor_effect[l] = linear_tumor_stimulus(
        // rep_vector(tumor_stim_pop_intercept + tumor_stim_location_intercept[l], n_grid_tumors), 
        rep_vector(0, n_grid_tumors), 
        rep_matrix(tumor_stim_pop_coef + tumor_stim_location_coef[l], n_grid_tumors), 
        tumor_covar[grid_tumors]
      );
    }
  }
  
  { // Some summary stats for survival
    tuple(vector[max_all_t], vector[max_all_t]) base_marginal_prob_res = calculate_marginal_dp_prob(base_pf_cond_prob, max_all_t);  
    tuple(vector[max_all_t], vector[max_all_t]) one_tumor_marginal_prob_res = calculate_marginal_dp_prob(one_tumor_pf_cond_prob, max_all_t);  
    
    base_dp_prob = base_marginal_prob_res.1';
    one_tumor_dp_prob = one_tumor_marginal_prob_res.1';
    
    base_survival = base_marginal_prob_res.2;
    one_tumor_survival = one_tumor_marginal_prob_res.2;
    
    base_dp_prob_not_censored = base_dp_prob / (1 - base_survival[max_all_t]);
    one_tumor_dp_prob_not_censored = one_tumor_dp_prob / (1 - one_tumor_survival[max_all_t]);
    
    base_cond_expected_pfs = base_dp_prob_not_censored * pfs_range_vec[:max_all_t];
    one_tumor_cond_expected_pfs = one_tumor_dp_prob_not_censored * pfs_range_vec[:max_all_t];
    
    base_cond_median_pfs = pfs_quantiles_from_prob(base_dp_prob_not_censored', { 0.5 })[1];
    one_tumor_cond_median_pfs = pfs_quantiles_from_prob(one_tumor_dp_prob_not_censored', { 0.5 })[1];
  }
}
