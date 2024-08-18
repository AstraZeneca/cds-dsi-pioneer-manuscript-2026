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
  #include "tumor/tumor_data.stan"
    
  array[n_patients] int<lower = 0> pfs; // How many periods after baseline did patient survive.
  array[n_patients] int<lower = 0> death_week; 
  array[n_patients] int<lower = 0, upper = 1> right_censored;

  // I use these to generate 2d hazard ratios over tumor sizes observed over time. 
  int<lower = 0> n_grid_tumors;
  array[n_grid_tumors] int<lower = 1, upper = sum(n_patient_tumors)> grid_tumors;
  
  // Hyperparam
  #include "baseline_hazard_hyperparam.stan"
  #include "tumor_stim_hyperparam.stan"
}

transformed data {
  #include "tumor/tumor_transformed_data.stan" 
  #include "pfs_transformed_data.stan"
  #include "tumor_stim_transformed_data.stan"
  
  if (min(n_patient_screening_t) <= 0) {
    reject("We need at least one screening measurement per patient, for now: ", min(n_patient_screening_t));
  }
}

parameters {
  // For joint modeling
  #include "tumor/tumor_parameters.stan"
  
  #include "baseline_hazard_parameters.stan"
  #include "tumor_stim_parameters.stan"
}

transformed parameters {
  #include "tumor/tumor_transformed_parameters.stan"
  #include "baseline_hazard_transformed_parameters.stan"
  #include "tumor_stim_transformed_parameters.stan"

  // This is a sum of the contribution of all a patient's tumors to their hazard 
  vector[n_patients] total_time_invar_tumor_stim;
  vector[n_patients] total_time_invar_tumor_stim_no_intercept; // This used outside the model, so don't delete it.
  vector[n_time_periods] disease_progress_pred; // DP predictor for all patients at all observed and unobserved intervals.
  
  (total_time_invar_tumor_stim, total_time_invar_tumor_stim_no_intercept, disease_progress_pred) = calc_disease_progress_pred_from_early_tumors(
    pfs, tumor_covar,
    tumor_hazard_type, n_patient_tumors,
    patient_trial,
    tumor_location,
    right_uncensored, interval_censored, 
    log_trial_lambda,
    tumor_stim_pop_intercept, tumor_stim_trial_intercept, tumor_stim_location_intercept, 
    tumor_stim_pop_coef, tumor_stim_trial_coef, tumor_stim_location_coef,
    gen_pfs, max_all_t
  ); 
  
  vector<lower = 0, upper = 1>[n_time_periods] disease_progress_prob = inv_cloglog(disease_progress_pred); // DP conditional probability for all patients (same as above)
}

model {
  // For joint modeling
  #include "tumor/tumor_model.stan"
  
  // Priors
  
  #include "baseline_hazard_priors.stan"
  #include "tumor_stim_priors.stan"

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
  
  if (tumor_hazard_type > 0) {
    one_tumor_pf_cond_prob = 
      1 - calculate_progress_linear_prob(log_lambda, tumor_stim_pop_intercept, tumor_stim_pop_coef, rep_row_vector(1, n_tumor_covar_col)); 
  } else {
    one_tumor_pf_cond_prob = base_pf_cond_prob;
  } 
  
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_base_pf_cond_prob; 
  array[add_trial_level ? n_trials : 0] vector<lower = 0, upper = 1>[max_all_t] trial_one_tumor_pf_cond_prob; 
  
  if (add_trial_level) {
    for (s in 1:n_trials) {
      trial_base_pf_cond_prob[s] = 1 - inv_cloglog(log_trial_lambda[s]);
      
      if (tumor_hazard_type > 0) {
        trial_one_tumor_pf_cond_prob[s] = 
          1 - calculate_progress_linear_prob(log_trial_lambda[s], tumor_stim_trial_intercept[s], tumor_stim_trial_coef[s], rep_row_vector(1, n_tumor_covar_col));
      } else {
        trial_one_tumor_pf_cond_prob[s] = trial_base_pf_cond_prob[s];
      } 
    }
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
