array[] int forecast_pfs_rng(
  int pred_calendar_week, array[] int start_calendar_week,
  array[] int patient_ids,
  array[] int pfs_interval_pos, array[] int pfs, array[] int pfs_right_censored, matrix log_cond_prob_surv,
  array[] int conf_resp_interval_pos, array[] int conf_resp_week, array[] int conf_resp, array[] int conf_resp_censored, matrix log_crcr_cond_prob_surv
) {
  int n_sample_patients = size(patient_ids);
  array[n_sample_patients] int pred_pfs = pfs[patient_ids];
  
  for (i in 1:n_sample_patients) {
    if (pfs_right_censored[patient_ids[i]] || (start_calendar_week[i] + pred_pfs[i] - 1) > pred_calendar_week) {
      int conf_resp_cause = conf_resp[patient_ids[i]] + 1;
      
      int max_weeks_observed = pred_calendar_week - start_calendar_week[i] + 1;
      
      if (conf_resp_censored[patient_ids[i]] || (start_calendar_week[i] + conf_resp_week[patient_ids[i]] - 1) > pred_calendar_week) {
        int n_weeks_obs_unclassified = min(conf_resp_week[patient_ids[i]], max_weeks_observed); 
        int curr_conf_resp_interval_pos = conf_resp_interval_pos[patient_ids[i]] + n_weeks_obs_unclassified;
        int curr_conf_resp_interval_end = conf_resp_interval_pos[patient_ids[i] + 1] - 1;
       
        conf_resp_cause = competing_risks_survival_time_rng(log_crcr_cond_prob_surv[curr_conf_resp_interval_pos:curr_conf_resp_interval_end]).3; 
      }
      
      int n_weeks_obs_surv = min(pred_pfs[i], max_weeks_observed); 
      int curr_pfs_interval_pos = pfs_interval_pos[patient_ids[i]] + n_weeks_obs_surv;
      int curr_pfs_interval_end = pfs_interval_pos[patient_ids[i] + 1] - 1;
      
      pred_pfs[i] = n_weeks_obs_surv + survival_time_rng(log_cond_prob_surv[curr_pfs_interval_pos:curr_pfs_interval_end, conf_resp_cause]).1; 
    }
  }
  
  return pred_pfs;
}

tuple(int, real) bootstrap_median_pfs_rng(
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
  
  bs_pfs[:n_sample_mature] = forecast_pfs_rng(
    pred_calendar_week,
    bs_sorted_start_calendar_week[:n_sample_mature],
    bs_sample_idx[mature_sorted_idx[:n_sample_mature]],
    pfs_interval_pos, pfs, pfs_right_censored, log_cond_prob_surv,
    conf_resp_interval_pos, conf_resp_week, conf_resp, conf_resp_censored, log_crcr_cond_prob_surv
  );
  
  bs_pfs[(n_sample_mature + 1):] = forecast_pfs_rng(
    pred_calendar_week,
    bs_sorted_start_calendar_week[(n_sample_mature + 1):][remaining_sort_idx[:n_sample_immature]],
    remaining_patients[:n_sample_immature],
    pfs_interval_pos, pfs, pfs_right_censored, log_cond_prob_surv,
    conf_resp_interval_pos, conf_resp_week, conf_resp, conf_resp_censored, log_crcr_cond_prob_surv
  );
  
  return(n_sample_immature, survival_median(bs_pfs, max_all_t).1);
}