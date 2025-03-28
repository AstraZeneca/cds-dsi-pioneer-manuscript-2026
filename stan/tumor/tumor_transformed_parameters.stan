vector[n_patients] patient_tumor_gp_intercept_effect;

for (s in 1:n_trials) {
  int patient_pos, patient_end;
  (patient_pos, patient_end) = get_pos(trial_patient_pos, s);
  
  patient_tumor_gp_intercept_effect[patient_pos:patient_end] = 
    raw_patient_tumor_gp_intercept_effect[patient_pos:patient_end] * patient_tumor_gp_intercept_sd[min(s, n_tumor_separate_trials)];
}

vector[separate_trial_tumor_gp ? 0 : n_pop_unique_visits] pop_obs_tumor_gp;
vector[sum(n_trial_unique_visits)] trial_obs_tumor_gp;

{
  int n_curr_pop_unique_visits = separate_trial_tumor_gp ? 0 : n_pop_unique_visits;
 
  if (!separate_trial_tumor_gp) { 
    pop_obs_tumor_gp = calc_gp_pred(all_measure_t[pop_unique_visits_idx], tumor_mean[1], pop_tumor_gp_alpha[1], pop_tumor_gp_rho[1], delta, pop_tumor_gp_eta);
  }
  
  for (s in 1:n_trials) {
    int actual_s = min(s, n_tumor_separate_trials);
    
    int trial_gp_pos, trial_gp_end;
    (trial_gp_pos, trial_gp_end) = get_pos(trial_unique_visits_pos, s);
   
    if (separate_trial_tumor_gp) {
      trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = calc_gp_pred(
        all_measure_t[get_int_sub_array(trial_unique_visits_idx, trial_unique_visits_pos, s)], 
        tumor_mean[s], pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], delta, get_sub_vector(pop_tumor_gp_eta, trial_unique_visits_pos, s) 
      );
    } else {
      array[n_trial_unique_visits[s]] int curr_trial2pop_idx = get_int_sub_array(trial2pop_unique_visit_idx, trial_unique_visits_pos, s);
      trial_obs_tumor_gp[trial_gp_pos:trial_gp_end] = pop_obs_tumor_gp[curr_trial2pop_idx];
    }
  }
}