array[n_bootstrap_cr_maturity_rates] int<lower = 1> bs_cr_prediction_calendar_week = rep_array(1000000, n_bootstrap_cr_maturity_rates);
array[n_bootstrap_cr_maturity_rates] int n_bs_sample_cr_classified = rep_array(0, n_bootstrap_cr_maturity_rates); 
array[n_bootstrap_cr_maturity_rates] int n_bs_sample_cr_unclassified = rep_array(0, n_bootstrap_cr_maturity_rates); 
vector<lower = 0, upper = 1>[n_bootstrap_cr_maturity_rates] bs_cr_orr = rep_vector(0, n_bootstrap_cr_maturity_rates);
vector<lower = 0, upper = 1>[n_bootstrap_cr_maturity_rates] bs_cr_conf_resp_censored_prop = rep_vector(0, n_bootstrap_cr_maturity_rates);
vector<lower = 0>[n_bootstrap_cr_maturity_rates] bs_cr_median_pfs = rep_vector(max_all_t, n_bootstrap_cr_maturity_rates);
vector<lower = 0, upper = 1>[n_bootstrap_cr_maturity_rates] bs_cr_pfs_right_censored_prop = rep_vector(0, n_bootstrap_cr_maturity_rates);
array[n_bootstrap_cr_maturity_rates] int n_bs_cr_conf_resp_predicted = rep_array(0, n_bootstrap_cr_maturity_rates);
array[n_bootstrap_cr_maturity_rates] int n_bs_cr_pfs_predicted = rep_array(0, n_bootstrap_cr_maturity_rates);

array[n_bootstrap_pfs_maturity_rates] int<lower = 1> bs_pfs_prediction_calendar_week = rep_array(1000000, n_bootstrap_pfs_maturity_rates);
array[n_bootstrap_pfs_maturity_rates] int n_bs_sample_pfs_progressed = rep_array(0, n_bootstrap_pfs_maturity_rates); 
array[n_bootstrap_pfs_maturity_rates] int n_bs_sample_pfs_surviving = rep_array(0, n_bootstrap_pfs_maturity_rates); 
vector<lower = 0, upper = 1>[n_bootstrap_pfs_maturity_rates] bs_pfs_orr = rep_vector(0, n_bootstrap_pfs_maturity_rates);
vector<lower = 0, upper = 1>[n_bootstrap_pfs_maturity_rates] bs_pfs_conf_resp_censored_prop = rep_vector(0, n_bootstrap_pfs_maturity_rates);
vector<lower = 0>[n_bootstrap_pfs_maturity_rates] bs_pfs_median_pfs = rep_vector(max_all_t, n_bootstrap_pfs_maturity_rates);
vector<lower = 0, upper = 1>[n_bootstrap_pfs_maturity_rates] bs_pfs_pfs_right_censored_prop = rep_vector(0, n_bootstrap_pfs_maturity_rates);
array[n_bootstrap_pfs_maturity_rates] int n_bs_pfs_conf_resp_predicted = rep_array(0, n_bootstrap_pfs_maturity_rates);
array[n_bootstrap_pfs_maturity_rates] int n_bs_pfs_pfs_predicted = rep_array(0, n_bootstrap_pfs_maturity_rates);

profile("leave-one-trial-out bootstrap") {
  if (leave_out_trial > 0 && n_bootstrap_sample > 0 && n_bootstrap_cr_maturity_rates > 0) {
    array[n_bootstrap_sample] int bs_sample_idx;
    array[n_bootstrap_sample] int bs_start_calendar_week; 
    array[n_bootstrap_sample] int bs_cr_mature_calendar_week; 
    array[n_bootstrap_sample] int bs_pfs_mature_calendar_week; 
    
    int patient_pos = trial_patient_pos[leave_out_trial];
    int patient_end = trial_patient_pos[leave_out_trial + 1] - 1;
    
    bs_sample_idx = discrete_range_rng(rep_array(patient_pos, n_bootstrap_sample), rep_array(patient_end, n_bootstrap_sample)); 
    
    for (bsi in 1:n_bootstrap_sample) {
      bs_cr_mature_calendar_week[bsi] = sorted_experiment_start_week[patient_pos + bsi - 1] + confirmed_response_week[bs_sample_idx[bsi]] - 1;
      bs_pfs_mature_calendar_week[bsi] = sorted_experiment_start_week[patient_pos + bsi - 1] + pfs[bs_sample_idx[bsi]] - right_censored[bs_sample_idx[bsi]];
    }
    
    array[n_bootstrap_sample] int cr_mature_sorted_idx = sort_indices_asc(bs_cr_mature_calendar_week);
    array[n_bootstrap_sample] int pfs_mature_sorted_idx = sort_indices_asc(bs_pfs_mature_calendar_week);
    
    bs_start_calendar_week = sorted_experiment_start_week[patient_pos:(patient_pos + n_bootstrap_sample - 1)];
    
    bs_cr_mature_calendar_week = bs_cr_mature_calendar_week[cr_mature_sorted_idx];
    bs_pfs_mature_calendar_week = bs_pfs_mature_calendar_week[pfs_mature_sorted_idx];
    
    (n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, bs_cr_prediction_calendar_week, 
     bs_cr_orr, bs_cr_conf_resp_censored_prop, n_bs_cr_conf_resp_predicted, 
     bs_cr_median_pfs, bs_cr_pfs_right_censored_prop, n_bs_cr_pfs_predicted) = 
      leave_out_trial_bootstrap_rng(
        bootstrap_cr_maturity_rates, bs_sample_idx,
        bs_start_calendar_week,
        cr_mature_sorted_idx, bs_cr_mature_calendar_week,
        max_all_t,
        patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
        patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
      );
    
    (n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, bs_pfs_prediction_calendar_week, 
     bs_pfs_orr, bs_pfs_conf_resp_censored_prop, n_bs_pfs_conf_resp_predicted, 
     bs_pfs_median_pfs, bs_pfs_pfs_right_censored_prop, n_bs_pfs_pfs_predicted) = 
      leave_out_trial_bootstrap_rng(
        bootstrap_pfs_maturity_rates, bs_sample_idx,
        bs_start_calendar_week,
        pfs_mature_sorted_idx, bs_pfs_mature_calendar_week,
        max_all_t,
        patient_pfs_interval_pos, pfs, right_censored, log_cond_prob_surv,
        patient_conf_resp_interval_pos, confirmed_response_week, confirmed_response, confirmed_response_censored, log_crcr_cond_prob_surv
      );
  }
}