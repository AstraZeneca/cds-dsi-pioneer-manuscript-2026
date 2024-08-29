array[n_bootstrap_cr_maturity_rates] int<lower = 1> bs_cr_prediction_calendar_week = rep_array(1000000, n_bootstrap_cr_maturity_rates);
array[n_bootstrap_pfs_maturity_rates] int<lower = 1> bs_pfs_prediction_calendar_week = rep_array(1000000, n_bootstrap_pfs_maturity_rates);
array[n_bootstrap_cr_maturity_rates] int n_bs_sample_cr_classified = rep_array(0, n_bootstrap_cr_maturity_rates); 
array[n_bootstrap_cr_maturity_rates] int n_bs_sample_cr_unclassified = rep_array(0, n_bootstrap_cr_maturity_rates); 
array[n_bootstrap_pfs_maturity_rates] int n_bs_sample_pfs_progressed = rep_array(0, n_bootstrap_pfs_maturity_rates); 
array[n_bootstrap_pfs_maturity_rates] int n_bs_sample_pfs_surviving = rep_array(0, n_bootstrap_pfs_maturity_rates); 
vector<lower = 0>[n_bootstrap_cr_maturity_rates] bs_cr_median_pfs = rep_vector(max_all_t, n_bootstrap_cr_maturity_rates);
vector<lower = 0>[n_bootstrap_pfs_maturity_rates] bs_pfs_median_pfs = rep_vector(max_all_t, n_bootstrap_pfs_maturity_rates);

{
  array[n_bootstrap_sample] int bs_sample_idx;
  array[n_bootstrap_sample] int bs_start_calendar_week; 
  array[n_bootstrap_sample] int bs_cr_mature_calendar_week; 
  array[n_bootstrap_sample] int bs_pfs_mature_calendar_week; 
  vector[n_bootstrap_sample] bs_cr_maturity_rate;
  vector[n_bootstrap_sample] bs_pfs_maturity_rate;
  
  profile("leave-one-trial-out bootstrap") {
    if (leave_out_trial > 0 && n_bootstrap_sample > 0 && n_bootstrap_cr_maturity_rates > 0) {
      int patient_pos = trial_patient_pos[leave_out_trial];
      int patient_end = trial_patient_pos[leave_out_trial + 1] - 1;
      
      bs_sample_idx = discrete_range_rng(rep_array(patient_pos, n_bootstrap_sample), rep_array(patient_end, n_bootstrap_sample)); 
      
      for (bsi in 1:n_bootstrap_sample) {
        bs_cr_mature_calendar_week[bsi] = sorted_experiment_start_week[patient_pos + bsi - 1] + confirmed_response_week[bs_sample_idx[bsi]] - 1;
        bs_pfs_mature_calendar_week[bsi] = sorted_experiment_start_week[patient_pos + bsi - 1] + pfs[bs_sample_idx[bsi]] - right_censored[bs_sample_idx[bsi]];
      }
      
      array[n_bootstrap_sample] int cr_mature_sorted_idx = sort_indices_asc(bs_cr_mature_calendar_week);
      array[n_bootstrap_sample] int pfs_mature_sorted_idx = sort_indices_asc(bs_pfs_mature_calendar_week);
      
      bs_cr_mature_calendar_week = bs_cr_mature_calendar_week[cr_mature_sorted_idx];
      bs_pfs_mature_calendar_week = bs_pfs_mature_calendar_week[pfs_mature_sorted_idx];
      bs_start_calendar_week = sorted_experiment_start_week[patient_pos:(patient_pos + n_bootstrap_sample - 1)]; // [cr_mature_sorted_idx];
      
      int cr_mature_rate_pos = 1;
      int pfs_mature_rate_pos = 1;
      
      for (bsi in 1:n_bootstrap_sample) {
        int cr_patient_id = bs_sample_idx[cr_mature_sorted_idx[bsi]]; 
        int pfs_patient_id = bs_sample_idx[pfs_mature_sorted_idx[bsi]]; 
        
        bs_cr_maturity_rate[bsi] = 1.0 * (1 - confirmed_response_censored[cr_patient_id]) / n_bootstrap_sample; 
        bs_pfs_maturity_rate[bsi] = 1.0 * (1 - right_censored[pfs_patient_id]) / n_bootstrap_sample; 
        
        if (bsi > 1) {
          bs_cr_maturity_rate[bsi] += bs_cr_maturity_rate[bsi - 1];
          bs_pfs_maturity_rate[bsi] += bs_pfs_maturity_rate[bsi - 1];
        }
        
        if (cr_mature_rate_pos <= n_bootstrap_cr_maturity_rates && bs_cr_maturity_rate[bsi] >= bootstrap_cr_maturity_rates[cr_mature_rate_pos]) {
          n_bs_sample_cr_classified[cr_mature_rate_pos] = bsi;
          bs_cr_prediction_calendar_week[cr_mature_rate_pos] = bs_cr_mature_calendar_week[bsi]; 
          
          (n_bs_sample_cr_unclassified[cr_mature_rate_pos], bs_cr_median_pfs[cr_mature_rate_pos]) =
            bootstrap_median_pfs_rng(
              n_bs_sample_cr_classified[cr_mature_rate_pos],
              bs_cr_prediction_calendar_week[cr_mature_rate_pos],
              bs_sample_idx,
              bs_start_calendar_week[cr_mature_sorted_idx],
              cr_mature_sorted_idx,
              max_all_t,
              patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
              patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
            ); 
          
          cr_mature_rate_pos += 1;
        }
        
        if (pfs_mature_rate_pos <= n_bootstrap_pfs_maturity_rates && bs_pfs_maturity_rate[bsi] >= bootstrap_pfs_maturity_rates[pfs_mature_rate_pos]) {
          n_bs_sample_pfs_progressed[pfs_mature_rate_pos] = bsi;
          bs_pfs_prediction_calendar_week[pfs_mature_rate_pos] = bs_pfs_mature_calendar_week[bsi]; 
          
          (n_bs_sample_pfs_surviving[pfs_mature_rate_pos], bs_pfs_median_pfs[pfs_mature_rate_pos]) =
            bootstrap_median_pfs_rng(
              n_bs_sample_pfs_progressed[pfs_mature_rate_pos],
              bs_pfs_prediction_calendar_week[pfs_mature_rate_pos],
              bs_sample_idx,
              bs_start_calendar_week[pfs_mature_sorted_idx],
              pfs_mature_sorted_idx,
              max_all_t,
              patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
              patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
            ); 
          
          pfs_mature_rate_pos += 1;
        }
      }
    }
  }
}