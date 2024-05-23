vector[fit_data ? n_patients : 0] log_lik; 

if (fit_data) {
  // Do not ignore censoring when calculating this
  log_lik = calc_pch_loglik(pfs, right_uncensored, interval_censored, 0, disease_progress_prob, gen_pfs ? max_all_t : 0, patient_max_2nd_tumor_t);
}

array[gen_pfs ? n_patients : 0] int<lower = 0> rep_pfs;
array[gen_pfs ? n_patients : 0] int<lower = 0, upper = 1> rep_right_censored;
array[gen_pfs ? n_patients : 0] int<lower = 0> rep_interval_censored;

// Kaplan-Meier survival probability, aggregated over generated patients' data.  
vector<lower = 0, upper = 1>[gen_pfs ? max(t_measure) + 1 : 0] km_est; 
array[add_trial_level && gen_pfs ? n_trials : 0] vector<lower = 0, upper = 1>[max(t_measure) + 1] trial_km_est; 

vector<lower = 0, upper = 1>[gen_pfs ? n_patients : 0] rep_pfs_6mon;
vector<lower = 0, upper = 1>[gen_pfs ? n_patients : 0] rep_pfs_9mon;

if (gen_pfs) { // Retrodiction, generating simulated data.
  int tumor_pos = 1;
  int pfs_interval_pos = 1;
  int t_pos = 1;
  
  for (i in 1:n_patients) {
    int pfs_interval_end = pfs_interval_pos + max_all_t - 1;
    int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
    int t_end = t_pos + n_measures[tumor_pos] - 1; // This is for just one tumor
    int t_all_end = t_pos + sum(n_measures[tumor_pos:tumor_end]) - 1; // For all the patient's tumors 
   
    // tuple(int, int, int, int) pfs_res = pfs_rng(
    //   disease_progress_prob[pfs_interval_pos:pfs_interval_end], 
    //   gen_interval_censored ? t_measure[(t_pos + n_screening_t[tumor_pos]):t_end] : pfs_range_int
    // );
    
    // rep_interval_censored[i] = pfs_res.1;
    // rep_right_censored[i] = pfs_res.2;
    // rep_pfs[i] = pfs_res.3;
    
    int rep_actual_pfs;
    
    (rep_interval_censored[i], rep_right_censored[i], rep_pfs[i], rep_actual_pfs) = pfs_rng(
      disease_progress_prob[pfs_interval_pos:pfs_interval_end], 
      gen_interval_censored ? t_measure[(t_pos + n_screening_t[tumor_pos]):t_end] : pfs_range_int
    );
    
    int max_t = max(t_measure[t_pos:t_all_end]);
    
    if (max_t >= 9 * 4 || (!rep_right_censored[i] && rep_actual_pfs <= 9 * 4)) {
      rep_pfs_6mon[i] = rep_actual_pfs >= 6 * 4; // What about interval censoring?
      rep_pfs_9mon[i] = rep_actual_pfs >= 9 * 4; // What about interval censoring?
    } else {
      int n_intervals = 9 * 4 - rep_actual_pfs;
      tuple(vector[n_intervals], vector[n_intervals]) marginal_prob_res = 
        calculate_marginal_dp_prob(disease_progress_prob[(pfs_interval_pos + rep_actual_pfs):(pfs_interval_pos + 9 * 4 - 1)], n_intervals);  
      
      rep_pfs_9mon[i] = 1 - marginal_prob_res.2[9 * 4 - rep_actual_pfs]; 
      
      if (max_t >= 6 * 4 || (!right_censored[i] && rep_actual_pfs <= 6 * 4)) {
        rep_pfs_6mon[i] = rep_actual_pfs >= 6 * 4; // What about interval censoring?
      } else {
        rep_pfs_6mon[i] = 1 - marginal_prob_res.2[6 * 4 - rep_actual_pfs]; 
      }
    } 
    
    pfs_interval_pos = pfs_interval_end + 1;
    t_pos += sum(n_measures[tumor_pos:tumor_end]);
    tumor_pos = tumor_end + 1;
  }
  
  km_est = estimate_kaplan_meier(rep_pfs, rep_right_censored, max_all_t).1; 
  
  if (add_trial_level) {
    int trial_patient_pos = 1;
    
    for (s in 1:n_trials) {
      int trial_patient_end = trial_patient_pos + n_trial_patients[s] - 1;
      
      trial_km_est[s] = estimate_kaplan_meier(rep_pfs[trial_patient_pos:trial_patient_end], rep_right_censored[trial_patient_pos:trial_patient_end], max_all_t).1; 
      
      trial_patient_pos = trial_patient_end + 1;
    }
  }
}
