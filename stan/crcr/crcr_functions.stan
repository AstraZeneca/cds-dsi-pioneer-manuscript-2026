/** Calculate the marginal probability of exit/event at every interval, separately for each competing risk. 
 *
 * @param log_crcr_cond_prob_surv The matrix of log conditional probabilities of survival
 * @return Matrix of log marginal probabilities
 */
matrix exit_log_marginal_prob(matrix log_crcr_cond_prob_surv) {
  int max_confresp_week = rows(log_crcr_cond_prob_surv);
  int n_causes = cols(log_crcr_cond_prob_surv);
  matrix[max_confresp_week, n_causes] log_marg_prob; 
  
  for (t in 1:max_confresp_week) {
    log_marg_prob[t] = sum(log_crcr_cond_prob_surv[1:(t - 1)]) + log1m_exp(log_crcr_cond_prob_surv[t]); 
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
array[] matrix exit_log_marginal_prob(int n_patients, matrix log_crcr_cond_prob_surv, int max_confresp_week) {
  int n_causes = cols(log_crcr_cond_prob_surv);
  array[n_patients] matrix[max_confresp_week, n_causes] log_marg_prob; 
  
  for (i in 1:n_patients) { 
    int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
    log_marg_prob[i] = exit_log_marginal_prob(log_crcr_cond_prob_surv[patient_prob_pos:(patient_prob_pos + max_confresp_week - 1)]);
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

/** Generate a survival profile for a single patient.
 *
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @return tuple(last interval before exit, is patient right censored, exit to which competing risk)
 */
tuple(int, int, int) competing_risks_survival_time_rng(matrix log_cond_prob_surv) {
  int n_intervals = rows(log_cond_prob_surv);
  int n_causes = cols(log_cond_prob_surv);
  
  matrix[n_intervals, n_causes] cond_prob_exit = 1 - exp(log_cond_prob_surv); 
 
  int survival_time = 0;
  int exit_cause = n_causes;
  
  for (t in 1:n_intervals) {
    array[n_causes] int exit_causes = bernoulli_rng(cond_prob_exit[t]);
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
tuple(int, int, int) competing_risks_survival_time_rng(matrix log_cond_prob_surv, int exit_cause, int event_time, int right_censored, int interval_censored) {
  int survival_time = event_time, forecast_right_censored = right_censored, forecast_exit_cause = exit_cause;
  
  if (right_censored) { // If right censored, forecast survival profile.
    (survival_time, forecast_right_censored, forecast_exit_cause) = competing_risks_survival_time_rng(log_cond_prob_surv[(event_time + 1):]);
    survival_time += event_time;
  } else if (interval_censored > 0) { // If interval censored, draw a survival time within the range of intervals.
    survival_time += interval_censored_survival_time_rng(log_cond_prob_surv[(event_time + 1):(event_time + interval_censored + 1), exit_cause]); 
  }
  
  return(survival_time, forecast_right_censored, forecast_exit_cause);
}

/** Calculate the survival log likelihood for all the patients.
 *
 * @param last_unclass_week Last week observed with no exit (survival).
 * @param event_cause What risk caused their exit.
 * @param right_censored
 * @param interval_censored
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param max_confresp_week The number of intervals for each patient 
 * @return Patient level log likelihoods
 */
vector calc_comp_risk_pch_loglik(
  array[] int last_unclass_week, array[] int event_cause, array[] int right_censored, array[] int interval_censored, matrix log_cond_prob_surv, int max_confresp_week
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
    
    real reuse_lp = sum(log_cond_prob_surv[interval_pos:interval_end]); // loglik for the known survival part
    
    // For IC, I need to create a mixture of all the possible true intervals of exit.
    vector[interval_censored[i] + 1] ic_mix_lp = rep_vector(reuse_lp, interval_censored[i] + 1);
    
    for (c in 0:interval_censored[i]) {
      ic_mix_lp[c + 1] += 
        sum(log_cond_prob_surv[(interval_end + 1):(interval_end + c)]) + 
        (1 - right_censored[i]) * log1m_exp(log_cond_prob_surv[interval_end + c + 1, event_cause[i]]);
    }
    
    lp[i] = log_sum_exp(ic_mix_lp) - log(interval_censored[i] + 1); 
    
    interval_pos += max_confresp_week; 
  }
  
  return lp;
}
 
/** Stan distribution _lpmf giving the sum of all patients' log likelihood. 
 *
 * @param last_unclass_week Last week observed with no exit (survival).
 * @param event_cause What risk caused their exit.
 * @param right_censored
 * @param interval_censored
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param max_confresp_week The number of intervals for each patient 
 * @return Patient level log likelihoods
 */
real comp_risk_pch_lpmf(
  array[] int last_unclass_week, array[] int event_cause, array[] int right_censored, array[] int interval_censored, matrix log_cond_prob_surv, int max_confresp_week
) {
  return sum(calc_comp_risk_pch_loglik(last_unclass_week, event_cause, right_censored, interval_censored, log_cond_prob_surv, max_confresp_week));
}

/** Calcuate the log likelihood split over multiple threads using reduce_sum().
 *
 * @param last_unclass_week Last week observed with no exit (survival).
 * @param start Index of first patient in this batch.
 * @param end Index of last patient in this batch.
 * @param event_cause What risk caused their exit.
 * @param right_censored
 * @param interval_censored
 * @param log_cond_prob_surv The matrix of log conditional probabilities of survival
 * @param max_confresp_week The number of intervals for each patient 
 * @return Patient level log likelihoods
 */
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
