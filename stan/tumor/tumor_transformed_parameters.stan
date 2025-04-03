vector[n_trials] trial_tumor_gp_intercept_effect = zeros_vector(n_trials);
vector[n_patients] patient_tumor_gp_intercept_effect;

vector[n_pop_unique_visits] pop_obs_tumor_gp;
vector[sum(n_trial_unique_visits)] trial_obs_tumor_gp;
vector[sum(n_patient_unique_visits)] patient_obs_tumor_gp;

if (separate_trial_tumor_gp) { 
  pop_obs_tumor_gp = zeros_vector(n_pop_unique_visits);
} else if (patient_gp_only) {
  pop_obs_tumor_gp = rep_vector(tumor_mean[1], n_pop_unique_visits);
} else {
  // pop_obs_tumor_gp = ncp_gp_matern52(
  pop_obs_tumor_gp = calc_gp_pred(
    all_tumor_measure_t[pop_unique_visits_idx], tumor_mean[1], pop_tumor_gp_alpha[1], pop_tumor_gp_rho[1], delta, pop_tumor_gp_eta
  );
} 

for (s in 1:n_trials) {
  int actual_s = min(s, n_tumor_separate_trials);
  
  int trial_gp_pos, trial_gp_end;
  (trial_gp_pos, trial_gp_end) = get_pos(trial_unique_visits_pos, s);
 
  if (separate_trial_tumor_gp) {
    trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = calc_gp_pred(
    // trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = ncp_gp_matern52(
      all_tumor_measure_t[get_int_sub_array(trial_unique_visits_idx, trial_unique_visits_pos, s)], 
      tumor_mean[s], pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], delta, get_sub_vector(pop_tumor_gp_eta, trial_unique_visits_pos, s) 
    );
  } else {
    array[n_trial_unique_visits[s]] int trial2pop_idx = get_int_sub_array(trial2pop_unique_visit_idx, trial_unique_visits_pos, s);
    
    if (add_trial_level_tumor) {
      trial_tumor_gp_intercept_effect[s] = raw_trial_tumor_gp_intercept_effect[s] * trial_tumor_gp_intercept_sd;
     
      if (patient_gp_only) { 
        trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = pop_obs_tumor_gp[trial2pop_idx] + trial_tumor_gp_intercept_effect[s];
      } else {
        trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = calc_gp_pred(
        // trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = ncp_gp_matern52(
          all_tumor_measure_t[get_int_sub_array(trial_unique_visits_idx, trial_unique_visits_pos, s)], 
          pop_obs_tumor_gp[trial2pop_idx] + trial_tumor_gp_intercept_effect[s], 
          trial_tumor_gp_alpha, trial_tumor_gp_rho, delta, get_sub_vector(trial_tumor_gp_eta, trial_unique_visits_pos, s) 
        );
      } 
    } else {
      trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = pop_obs_tumor_gp[trial2pop_idx];
    }
  } 
  
  int patient_pos, patient_end;
  (patient_pos, patient_end) = get_pos(trial_patient_pos, s);
  
  patient_tumor_gp_intercept_effect[patient_pos:patient_end] = raw_patient_tumor_gp_intercept_effect[patient_pos:patient_end] * patient_tumor_gp_intercept_sd[actual_s];
  
  for (i in patient_pos:patient_end) {
    int patient_gp_pos, patient_gp_end;
    (patient_gp_pos, patient_gp_end) = get_pos(patient_unique_visits_pos, i);
    
    array[n_patient_unique_visits[i]] int patient2trial_idx = get_int_sub_array(patient2trial_unique_visit_idx, patient_unique_visits_pos, i);
    
    patient_obs_tumor_gp[patient_gp_pos:patient_gp_end] = calc_gp_pred(
    // patient_obs_tumor_gp[patient_gp_pos:patient_gp_end] = ncp_gp_matern52(
      all_tumor_measure_t[get_int_sub_array(patient_unique_visits_idx, patient_unique_visits_pos, i)], 
      trial_obs_tumor_gp[patient2trial_idx] + patient_tumor_gp_intercept_effect[i], 
      patient_tumor_gp_alpha[actual_s], patient_tumor_gp_rho[actual_s], delta, get_sub_vector(patient_tumor_gp_eta, patient_unique_visits_pos, i) 
    );
  }
}