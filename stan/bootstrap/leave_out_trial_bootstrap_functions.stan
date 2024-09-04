tuple(array[] int, array[] int, int, array[] int, array[] int, int) forecast_conf_resp_pfs_rng(
  int pred_calendar_week, array[] int start_calendar_week,
  array[] int patient_ids,
  array[] int pfs_interval_pos, array[] int pfs, array[] int pfs_right_censored, matrix log_cond_prob_surv,
  array[] int conf_resp_interval_pos, array[] int conf_resp_week, array[] int conf_resp, array[] int conf_resp_censored, matrix log_crcr_cond_prob_surv
) {
  int n_sample_patients = size(patient_ids);
  array[n_sample_patients] int pred_pfs = pfs[patient_ids];
  array[n_sample_patients] int pred_pfs_right_censored = pfs_right_censored[patient_ids];
  array[n_sample_patients] int pred_conf_resp_week = conf_resp_week[patient_ids];
  array[n_sample_patients] int pred_conf_resp = conf_resp[patient_ids];
  array[n_sample_patients] int pred_conf_resp_censored = conf_resp_censored[patient_ids];
  
  int n_conf_resp_predicted = 0, n_pfs_predicted = 0;
  
  for (i in 1:n_sample_patients) {
    int max_weeks_observed = pred_calendar_week - start_calendar_week[i] + 1;
    int conf_resp_cause = conf_resp[patient_ids[i]] + 1;
    
    if (pred_conf_resp_censored[i] || (start_calendar_week[i] + pred_conf_resp_week[i] - 1) > pred_calendar_week) {
      int n_weeks_obs_unclassified = min(pred_conf_resp_week[i], max_weeks_observed); 
      int curr_conf_resp_interval_pos = conf_resp_interval_pos[patient_ids[i]] + n_weeks_obs_unclassified;
      int curr_conf_resp_interval_end = conf_resp_interval_pos[patient_ids[i] + 1] - 1;
      
      int continued_survival;
     
      (continued_survival, pred_conf_resp_censored[i], conf_resp_cause) = 
        competing_risks_survival_time_rng(log_crcr_cond_prob_surv[curr_conf_resp_interval_pos:curr_conf_resp_interval_end]);
        
      pred_conf_resp_week[i] = n_weeks_obs_unclassified + continued_survival + (1 - pred_conf_resp_censored[i]);
      pred_conf_resp[i] = conf_resp_cause - 1;
      n_conf_resp_predicted += 1;
    }
    
    if (pred_pfs_right_censored[i] || (start_calendar_week[i] + pred_pfs[i] - 1) > pred_calendar_week) {
      int n_weeks_obs_surv = min(pred_pfs[i], max_weeks_observed); 
      
      int curr_pfs_interval_pos = pfs_interval_pos[patient_ids[i]] + n_weeks_obs_surv;
      int curr_pfs_interval_end = pfs_interval_pos[patient_ids[i] + 1] - 1;
      
      int continued_pfs;
      
      (continued_pfs, pred_pfs_right_censored[i]) = survival_time_rng(log_cond_prob_surv[curr_pfs_interval_pos:curr_pfs_interval_end, conf_resp_cause]); 
      
      pred_pfs[i] = n_weeks_obs_surv + continued_pfs; 
      n_pfs_predicted += 1;
    }
  }
  
  return(pred_conf_resp, pred_conf_resp_censored, n_conf_resp_predicted, pred_pfs, pred_pfs_right_censored, n_pfs_predicted);
}

tuple(int, real, real, int, real, real, int) bootstrap_orr_median_pfs_rng(
  int n_sample_mature, int pred_calendar_week, array[] int bs_sample_idx, array[] int bs_sorted_start_calendar_week, array[] int mature_sorted_idx,
  int max_all_t,
  array[] int pfs_interval_pos, array[] int pfs, array[] int pfs_right_censored, matrix log_cond_prob_surv,
  array[] int conf_resp_interval_pos, array[] int conf_resp_week, array[] int conf_resp, array[] int conf_resp_censored, matrix log_crcr_cond_prob_surv
) {
  int n_bootstrap_sample = size(bs_sample_idx);
  
  int n_remaining = n_bootstrap_sample - n_sample_mature;
  array[n_remaining] int remaining_patients, remaining_sort_idx;
  
  if (n_remaining > 0) {
    remaining_sort_idx = sort_indices_asc(bs_sorted_start_calendar_week[(n_sample_mature + 1):]);
    remaining_patients = bs_sample_idx[mature_sorted_idx[(n_sample_mature + 1):]][remaining_sort_idx];
  }
  
  int n_sample_immature = 0;
  
  for (ri in 1:n_remaining) {
    if (bs_sorted_start_calendar_week[(n_sample_mature + 1):][remaining_sort_idx[ri]] <= pred_calendar_week) {
      n_sample_immature += 1; 
    } else {
      break;
    }
  }
  
  int n_bs_sample = n_sample_mature + n_sample_immature;
  array[n_bs_sample] int bs_pfs;
  array[n_bs_sample] int bs_pfs_right_censored;
  int n_bs_pfs_predicted_mature, n_bs_pfs_predicted_immature;
  array[n_bs_sample] int bs_conf_resp;
  array[n_bs_sample] int bs_conf_resp_censored;
  int n_bs_conf_resp_predicted_mature, n_bs_conf_resp_predicted_immature;
  
  (bs_conf_resp[:n_sample_mature], bs_conf_resp_censored[:n_sample_mature], n_bs_conf_resp_predicted_mature, 
   bs_pfs[:n_sample_mature], bs_pfs_right_censored[:n_sample_mature], n_bs_pfs_predicted_mature) = forecast_conf_resp_pfs_rng(
    pred_calendar_week,
    bs_sorted_start_calendar_week[:n_sample_mature],
    bs_sample_idx[mature_sorted_idx[:n_sample_mature]],
    pfs_interval_pos, pfs, pfs_right_censored, log_cond_prob_surv,
    conf_resp_interval_pos, conf_resp_week, conf_resp, conf_resp_censored, log_crcr_cond_prob_surv
  );
  
  (bs_conf_resp[(n_sample_mature + 1):], bs_conf_resp_censored[(n_sample_mature + 1):], n_bs_conf_resp_predicted_immature, 
   bs_pfs[(n_sample_mature + 1):], bs_pfs_right_censored[(n_sample_mature + 1):], n_bs_pfs_predicted_immature) = forecast_conf_resp_pfs_rng(
    pred_calendar_week,
    bs_sorted_start_calendar_week[(n_sample_mature + 1):][remaining_sort_idx[:n_sample_immature]],
    remaining_patients[:n_sample_immature],
    pfs_interval_pos, pfs, pfs_right_censored, log_cond_prob_surv,
    conf_resp_interval_pos, conf_resp_week, conf_resp, conf_resp_censored, log_crcr_cond_prob_surv
  );
  
  return(n_sample_immature, 
         mean(bs_conf_resp), mean(bs_conf_resp_censored), n_bs_conf_resp_predicted_mature + n_bs_conf_resp_predicted_immature, 
         survival_median(bs_pfs, max_all_t).1, mean(bs_pfs_right_censored), n_bs_pfs_predicted_mature + n_bs_pfs_predicted_immature);
}

tuple(array[] int, array[] int, array[] int, vector, vector, array[] int, vector, vector, array[] int) leave_out_trial_bootstrap_rng(
  vector bootstrap_maturity_rates, array[] int bs_sample_idx,
  array[] int bs_start_calendar_week,
  array[] int mature_sorted_idx, array[] int bs_mature_calendar_week, 
  int max_all_t,
  array[] int pfs_interval_pos, array[] int pfs, array[] int pfs_right_censored, matrix log_cond_prob_surv,
  array[] int conf_resp_interval_pos, array[] int conf_resp_week, array[] int conf_resp, array[] int conf_resp_censored, matrix log_crcr_cond_prob_surv
) {
  int n_bootstrap_sample = size(bs_sample_idx);
  int n_bootstrap_maturity_rates = size(bootstrap_maturity_rates); 
  
  vector[n_bootstrap_sample] bs_maturity_rate;
  array[n_bootstrap_maturity_rates] int n_bs_sample_observed = rep_array(0, n_bootstrap_maturity_rates); 
  array[n_bootstrap_maturity_rates] int n_bs_sample_unobserved = rep_array(0, n_bootstrap_maturity_rates);
  array[n_bootstrap_maturity_rates] int bs_prediction_calendar_week = rep_array(1000000, n_bootstrap_maturity_rates); 
  vector[n_bootstrap_maturity_rates] bs_orr = rep_vector(0, n_bootstrap_maturity_rates);
  vector[n_bootstrap_maturity_rates] bs_conf_resp_censored_prop = rep_vector(1, n_bootstrap_maturity_rates);
  array[n_bootstrap_maturity_rates] int n_bs_conf_resp_predicted = rep_array(0, n_bootstrap_maturity_rates);
  vector[n_bootstrap_maturity_rates] bs_median_pfs = rep_vector(max_all_t, n_bootstrap_maturity_rates);
  vector[n_bootstrap_maturity_rates] bs_pfs_right_censored_prop = rep_vector(1, n_bootstrap_maturity_rates);
  array[n_bootstrap_maturity_rates] int n_bs_pfs_predicted = rep_array(0, n_bootstrap_maturity_rates);
  
  int mature_rate_pos = 1;
  int bsi = 1;
  
  while (bsi <= n_bootstrap_sample && mature_rate_pos <= n_bootstrap_maturity_rates) {
    int patient_id = bs_sample_idx[mature_sorted_idx[bsi]]; 
    
    bs_maturity_rate[bsi] = 1.0 * (1 - conf_resp_censored[patient_id]) / n_bootstrap_sample; 
    
    if (bsi > 1) {
      bs_maturity_rate[bsi] += bs_maturity_rate[bsi - 1];
    }
    
    int already_calculated = 0;
    
    while (mature_rate_pos <= n_bootstrap_maturity_rates && bs_maturity_rate[bsi] >= bootstrap_maturity_rates[mature_rate_pos]) {
      n_bs_sample_observed[mature_rate_pos] = bsi;
      bs_prediction_calendar_week[mature_rate_pos] = bs_mature_calendar_week[bsi]; 
      
      if (!already_calculated) { 
        (n_bs_sample_unobserved[mature_rate_pos], 
         bs_orr[mature_rate_pos], bs_conf_resp_censored_prop[mature_rate_pos], n_bs_conf_resp_predicted[mature_rate_pos], 
         bs_median_pfs[mature_rate_pos], bs_pfs_right_censored_prop[mature_rate_pos], n_bs_pfs_predicted[mature_rate_pos]) =
          bootstrap_orr_median_pfs_rng(
            n_bs_sample_observed[mature_rate_pos],
            bs_prediction_calendar_week[mature_rate_pos],
            bs_sample_idx,
            bs_start_calendar_week[mature_sorted_idx],
            mature_sorted_idx,
            max_all_t,
            pfs_interval_pos, pfs, pfs_right_censored, log_cond_prob_surv,
            conf_resp_interval_pos, conf_resp_week, conf_resp, conf_resp_censored, log_crcr_cond_prob_surv
          ); 
          
          already_calculated = 1;
      } else {
        (n_bs_sample_unobserved[mature_rate_pos], 
         bs_orr[mature_rate_pos], bs_conf_resp_censored_prop[mature_rate_pos], n_bs_conf_resp_predicted[mature_rate_pos], 
         bs_median_pfs[mature_rate_pos], bs_pfs_right_censored_prop[mature_rate_pos], n_bs_pfs_predicted[mature_rate_pos]) =
           (n_bs_sample_unobserved[mature_rate_pos - 1], 
           bs_orr[mature_rate_pos - 1], bs_conf_resp_censored_prop[mature_rate_pos - 1], n_bs_conf_resp_predicted[mature_rate_pos - 1], 
           bs_median_pfs[mature_rate_pos - 1], bs_pfs_right_censored_prop[mature_rate_pos - 1], n_bs_pfs_predicted[mature_rate_pos - 1]);
      }
      
      mature_rate_pos += 1;
    }
    
    bsi += 1;
  }
  
  return(
    n_bs_sample_observed, n_bs_sample_unobserved, bs_prediction_calendar_week, 
    bs_orr, bs_conf_resp_censored_prop, n_bs_conf_resp_predicted, 
    bs_median_pfs, bs_pfs_right_censored_prop, n_bs_pfs_predicted
  );
}