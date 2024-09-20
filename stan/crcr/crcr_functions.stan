matrix exit_log_marginal_prob(matrix log_crcr_cond_prob_surv) {
  int max_confresp_week = rows(log_crcr_cond_prob_surv);
  int n_causes = cols(log_crcr_cond_prob_surv);
  matrix[max_confresp_week, n_causes] log_marg_prob; 
  
  for (t in 1:max_confresp_week) {
    log_marg_prob[t] = sum(log_crcr_cond_prob_surv[1:(t - 1)]) + log1m_exp(log_crcr_cond_prob_surv[t]); 
  }
  
  return(log_marg_prob);
}

array[] matrix exit_log_marginal_prob(int n_patients, matrix log_crcr_cond_prob_surv, int max_confresp_week) {
  int n_causes = cols(log_crcr_cond_prob_surv);
  array[n_patients] matrix[max_confresp_week, n_causes] log_marg_prob; 
  
  for (i in 1:n_patients) { 
    int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
    log_marg_prob[i] = exit_log_marginal_prob(log_crcr_cond_prob_surv[patient_prob_pos:(patient_prob_pos + max_confresp_week - 1)]);
  }
  
  return(log_marg_prob);
}

tuple(array[] matrix, matrix) calc_cif(int n_patients, matrix log_crcr_cond_prob_surv, int max_confresp_week) {
  int n_causes = cols(log_crcr_cond_prob_surv);
  
  array[n_patients] matrix[max_confresp_week, n_causes] log_marg_prob = exit_log_marginal_prob(n_patients, log_crcr_cond_prob_surv, max_confresp_week); 
  
  array[n_patients] matrix[max_confresp_week, n_causes] cif; // cumulative incidence function
  matrix[n_patients, n_causes] prob_cause; 
  
  for (i in 1:n_patients) { 
    for (k in 1:n_causes) {
      cif[i, , k] = cumulative_sum(exp(log_marg_prob[i, , k]));
    }
    
    prob_cause[i] = cif[i, max_confresp_week];
    
    prob_cause[i] /= sum(prob_cause[i]); 
  }
  
  return(cif, prob_cause);
}

tuple(int, int, int) competing_risks_survival_time_rng(matrix log_cond_prob_surv) {
  int n_intervals = rows(log_cond_prob_surv);
  int n_causes = cols(log_cond_prob_surv);
  
  matrix[n_intervals, n_causes] cond_prob_exit = 1 - exp(log_cond_prob_surv); 
 
  int survival_time = 0;
  int exit_cause = n_causes;
  
  for (t in 1:n_intervals) {
    array[n_causes] int exit_causes = bernoulli_rng(cond_prob_exit[t]);
    int num_exits = sum(exit_causes);
    
    // exit_cause = categorical_rng(cond_prob_exit[t]');
  
    if (num_exits == 0) {
    // if (exit_cause > n_causes) {
      survival_time += 1;
    } else {
      if (num_exits == 1) {
        exit_cause = sort_indices_desc(exit_causes)[1];
      } else {
        exit_cause = sort_indices_desc(exit_causes)[discrete_range_rng(1, num_exits)];
      }
      
      break;
    }
  }
  
  return(survival_time, survival_time >= n_intervals, exit_cause);
}

vector calc_comp_risk_pch_loglik(
  array[] int last_unclass_week,
  array[] int event_cause,
  array[] int right_censored,
  array[] int interval_censored,
  matrix log_cond_prob_surv,
  data int max_confresp_week
) 
{
  int n_patients = size(last_unclass_week);
  int n_causes = cols(log_cond_prob_surv);
  vector[n_patients] lp;
  
  int interval_pos = 1;
    
  for (i in 1:n_patients) {
    if (right_censored[i] && interval_censored[i] > 0) {
      fatal_error("Interval censoring not allowed with right censored observations.");
    }
    
    int interval_end = interval_pos + last_unclass_week[i] - 1; 
    
    real reuse_lp = sum(log_cond_prob_surv[interval_pos:interval_end]);
    
    vector[interval_censored[i] + 1] ic_mix_lp = rep_vector(reuse_lp - log(interval_censored[i] + 1), interval_censored[i] + 1);
    
    for (c in 0:interval_censored[i]) {
      ic_mix_lp[c + 1] += 
        sum(log_cond_prob_surv[(interval_end + 1):(interval_end + c)]) + 
        (1 - right_censored[i]) * log1m_exp(log_cond_prob_surv[interval_end + c + 1, event_cause[i]]);
    }
    
    lp[i] = log_sum_exp(ic_mix_lp); 
    
    interval_pos += max_confresp_week; 
  }
  
  return lp;
}
 
real comp_risk_pch_lpmf(
  array[] int last_unclass_week, array[] int event_cause, array[] int right_censored, array[] int interval_censored, matrix log_cond_prob_surv, int max_confresp_week
) {
  return sum(calc_comp_risk_pch_loglik(last_unclass_week, event_cause, right_censored, interval_censored, log_cond_prob_surv, max_confresp_week));
}

real partial_sum_crcr_lpmf(
  array[] int last_unclass_week, int start, int end, 
  array[] int event_cause, array[] int right_censored, array[] int interval_censored, matrix log_cond_prob_surv, int max_confresp_week
) {
  int patient_interval_pos = 1 + (start - 1) * max_confresp_week; 
  int patient_interval_end = end * max_confresp_week; 
  
  return(comp_risk_pch_lpmf(
    last_unclass_week | event_cause[start:end], 
                        right_censored[start:end], interval_censored[start:end], 
                        log_cond_prob_surv[patient_interval_pos:patient_interval_end], max_confresp_week
  ));
}
