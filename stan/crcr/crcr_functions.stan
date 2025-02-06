/** Calculate the marginal probability of exit/event at every interval, separately for each competing risk. 
 *
 * @param log_crcr_cond_prob_surv The matrix of log conditional probabilities of survival
 * @return Matrix of log marginal probabilities
 */
array[] row_vector exit_log_marginal_prob(array[] row_vector log_crcr_cond_prob_surv) {
  int n_causes = size(log_crcr_cond_prob_surv), max_confresp_week = cols(log_crcr_cond_prob_surv[1]);
  
  array[n_causes] row_vector[max_confresp_week] log_marg_prob; 
  
  for (k in 1:n_causes) {
    for (t in 1:max_confresp_week) {
      log_marg_prob[k, t] = sum(log_crcr_cond_prob_surv[k, 1:(t - 1)]) + log1m_exp(log_crcr_cond_prob_surv[k, t]); 
    }
  }
  
  return(log_marg_prob);
}

/** Calculate the marginal probability of exit/event at every interval, separately for each competing risk, for all patients.  
 *
 * @param n_patients Number of patients
 * @param log_crcr_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param max_confresp_week The number of intervals for each patient 
 * @return Array of matrices of log marginal probabilities
 */
array[] matrix exit_log_marginal_prob(array[] matrix log_crcr_cond_prob_surv) {
  int n_causes = size(log_crcr_cond_prob_surv), n_patients = rows(log_crcr_cond_prob_surv[1]), max_confresp_week = cols(log_crcr_cond_prob_surv[1]);
  array[n_causes] matrix[n_patients, max_confresp_week] log_marg_prob; 
  
  for (i in 1:n_patients) { 
    log_marg_prob[, i] = exit_log_marginal_prob(log_crcr_cond_prob_surv[, i]);
  }
  
  return(log_marg_prob);
}

/** Calculate the cumulative incidence function (CIF) per patient and competing risk.
 *
 * @param n_patients Number of patients
 * @param log_crcr_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param max_confresp_week The number of intervals for each patient 
 * @return tuple of patient level CIF information and the approximate probability of each patient exiting to each of the competing risks. 
 */
tuple(array[] matrix, matrix) calc_cif(array[] matrix log_crcr_cond_prob_surv) {
  int n_causes = size(log_crcr_cond_prob_surv), n_patients = rows(log_crcr_cond_prob_surv[1]), max_confresp_week = cols(log_crcr_cond_prob_surv[1]);
  
  array[n_causes] matrix[n_patients, max_confresp_week] log_marg_prob = exit_log_marginal_prob(log_crcr_cond_prob_surv); 
  
  array[n_causes] matrix[n_patients, max_confresp_week] cif; // cumulative incidence function
  matrix[n_causes, n_patients] prob_cause; 
  
  for (i in 1:n_patients) { 
    for (k in 1:n_causes) {
      cif[k, i] = cumulative_sum(exp(log_marg_prob[k, i]));
    
      prob_cause[k, i] = cif[k, i, max_confresp_week];
    }
    
    prob_cause[, i] /= sum(prob_cause[, i]); 
  }
  
  return(cif, prob_cause);
}

/** Generate a survival profile for a single patient.
 *
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @return tuple(last interval before exit, is patient right censored, exit to which competing risk)
 */
tuple(int, int, int) competing_risks_survival_time_rng(array[] row_vector log_cond_prob_surv) {
  int n_intervals = cols(log_cond_prob_surv[1]);
  int n_causes = size(log_cond_prob_surv);
  
  matrix[n_causes, n_intervals] cond_prob_exit = 1 - exp(to_matrix(log_cond_prob_surv));
 
  int survival_time = 0;
  int exit_cause = n_causes;
  
  for (t in 1:n_intervals) {
    array[n_causes] int exit_causes = bernoulli_rng(cond_prob_exit[, t]);
    int num_exits = sum(exit_causes);
  
    if (num_exits == 0) { // Didn't exit to any of the competing risks
      survival_time += 1;
    } else {
      if (num_exits == 1) {
        exit_cause = sort_indices_desc(exit_causes)[1];
      } else { // Multiple candidate risks: randomly pick one.
        exit_cause = sort_indices_desc(exit_causes)[discrete_range_rng(1, num_exits)];
      }
      
      break;
    }
  }
  
  return(survival_time, survival_time >= n_intervals, exit_cause);
}

/** Generate a forecast survival profile for a single patient.
 *
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param exit_cause Observed exit cause
 * @param event_time When was exit observed
 * @param right_censored
 * @param interval_censored
 * @return tuple(last interval before exit, is patient right censored, exit to which competing risk)
 */
tuple(int, int, int) competing_risks_survival_time_rng(array[] row_vector log_cond_prob_surv, int exit_cause, int event_time, int right_censored, int interval_censored) {
  int survival_time = event_time, forecast_right_censored = right_censored, forecast_exit_cause = exit_cause;
  
  if (right_censored) { // If right censored, forecast survival profile.
    (survival_time, forecast_right_censored, forecast_exit_cause) = competing_risks_survival_time_rng(log_cond_prob_surv[, (event_time + 1):]);
    survival_time += event_time;
  } else if (interval_censored > 0) { // If interval censored, draw a survival time within the range of intervals.
    survival_time += interval_censored_survival_time_rng(log_cond_prob_surv[exit_cause, (event_time + 1):(event_time + interval_censored + 1)]); 
  }
  
  return(survival_time, forecast_right_censored, forecast_exit_cause);
}

real partial_sum_crcr_lpmf(
  array[] int last_unclass_week, int start, int end, 
  array[] int event_cause, array[] int right_censored, array[] int interval_censored, array[] matrix log_cond_prob_surv
) {
  return pch_lpmf(last_unclass_week | event_cause[start:end], right_censored[start:end], interval_censored[start:end], 0, log_cond_prob_surv[, start:end]);
}

