/**
 * Calculate Log Marginal Probabilities of Exit for Competing Risks
 *
 * This function computes the log marginal probabilities of exit for each cause
 * in a competing risks scenario, based on the conditional survival probabilities.
 *
 * The function performs the following steps:
 * 1. Calculates the overall survival probability for each patient at each time point.
 * 2. Computes the marginal probability of exit for each cause, patient, and time point.
 *
 * The marginal probability of exit for cause k at time t is calculated as:
 * P(Exit due to cause k at time t) = P(Survive until t-1) * P(Exit due to cause k at t | Survived until t-1)
 *
 * Note: All probabilities are computed and returned in log scale for numerical stability.
 *
 * @param log_crcr_cond_prob_surv Array of matrices containing log conditional 
 *        survival probabilities for each cause, patient, and time point.
 *        Dimensions: [n_causes, n_patients, max_confresp_week]
 *
 * @return Array of matrices containing log marginal probabilities of exit
 *         for each cause, patient, and time point.
 *         Dimensions: [n_causes, n_patients, max_confresp_week]
 */
array[] matrix exit_log_marginal_prob(array[] matrix log_crcr_cond_prob_surv) {
  int n_causes = size(log_crcr_cond_prob_surv);
  int n_patients = rows(log_crcr_cond_prob_surv[1]);
  int max_confresp_week = cols(log_crcr_cond_prob_surv[1]);
  
  array[n_causes] matrix[n_patients, max_confresp_week] log_marg_prob;
  
  for (i in 1:n_patients) {
    vector[max_confresp_week] log_surv = zeros_vector(max_confresp_week);
  
    // Calculate overall survival probability
    for (t in 1:max_confresp_week) {
      log_surv[t] = sum(log_crcr_cond_prob_surv[, i, t]);
  
      if (t > 1) {  
        log_surv[t] += log_surv[t-1];
      }
    }
  
    // Calculate marginal probability for each cause
    for (k in 1:n_causes) {
      for (t in 1:max_confresp_week) {
        log_marg_prob[k, i, t] = log1m_exp(log_crcr_cond_prob_surv[k, i, t]);
  
        if (t > 1) {  
          log_marg_prob[k, i, t] += log_surv[t-1];
        }
      }
    }
  }
  
  return log_marg_prob;
}

/**
 * Calculate the Cumulative Incidence Function (CIF) for Competing Risks
 *
 * This function computes the log of the Cumulative Incidence Function (CIF) for each patient
 * and competing risk, based on the conditional survival probabilities.
 *
 * @param log_crcr_cond_prob_surv Array of matrices containing log conditional 
 *        survival probabilities for each cause, patient, and time point.
 *        Dimensions: [n_causes, n_patients, max_confresp_week]
 *
 * @return Array of matrices containing the log of the Cumulative Incidence Function
 *         for each cause, patient, and time point.
 *         Dimensions: [n_causes, n_patients, max_confresp_week]
 *
 * The function performs the following steps:
 * 1. Calculates the log marginal probabilities of exit using the exit_log_marginal_prob function.
 * 2. Computes the log CIF for each cause, patient, and time point.
 *
 * The CIF for cause k at time t is calculated as:
 * CIF_k(t) = sum_{s=1}^t P(Exit due to cause k at time s)
 *
 * In log scale, this becomes:
 * log(CIF_k(t)) = log_sum_exp(log(P(Exit due to cause k at time s))) for s = 1 to t
 *
 * Note: All probabilities are computed and returned in log scale for numerical stability.
 */
array[] matrix calc_log_cif(array[] matrix log_crcr_cond_prob_surv) {
  int n_causes = size(log_crcr_cond_prob_surv);
  int n_patients = rows(log_crcr_cond_prob_surv[1]);
  int max_confresp_week = cols(log_crcr_cond_prob_surv[1]);
  
  array[n_causes] matrix[n_patients, max_confresp_week] log_marg_prob = exit_log_marginal_prob(log_crcr_cond_prob_surv);  
  
  array[n_causes] matrix[n_patients, max_confresp_week] log_cif;
  
  for (i in 1:n_patients) {
    for (k in 1:n_causes) {
      for (t in 1:max_confresp_week) {
        log_cif[k, i, t] = log_sum_exp(log_marg_prob[k, i, :t]);
      }
    }
  }
  
  return log_cif;
}

/**
 * Generate a survival profile for a single patient with competing risks
 *
 * This function simulates a survival time and exit cause for a patient in a 
 * competing risks scenario, based on the provided conditional survival probabilities.
 *
 * @param log_cond_prob_surv Array of row vectors containing log conditional 
 *        probabilities of survival for each cause and time interval.
 *        Dimensions: [n_causes, n_intervals]
 *
 * @return A tuple containing:
 *         1. int: The last interval before exit (survival time)
 *         2. int: Indicator if the patient is right censored (1 if censored, 0 otherwise)
 *         3. int: The cause of exit (1 to n_causes)
 *
 * The function simulates the survival process by:
 * 1. Converting conditional survival probabilities to exit probabilities
 * 2. For each time interval, simulating potential exits for all causes
 * 3. If multiple exits occur in the same interval, randomly selecting one
 * 4. Continuing until an exit occurs or all intervals are exhausted
 */
tuple(int, int, int) competing_risks_survival_time_rng(array[] row_vector log_cond_prob_surv) {
  int n_intervals = cols(log_cond_prob_surv[1]);
  int n_causes = size(log_cond_prob_surv);
  
  matrix[n_causes, n_intervals] mat_log_cond_prov_surv = to_matrix(log_cond_prob_surv);
  matrix[n_causes, n_intervals] log_odds_cond_prob_exit = log1m_exp(mat_log_cond_prov_surv) - mat_log_cond_prov_surv;
 
  int survival_time = 0;
  int exit_cause = n_causes;
  
  for (t in 1:n_intervals) {
    array[n_causes] int exit_causes = bernoulli_logit_rng(log_odds_cond_prob_exit[, t]);
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

/**
 * Generate a forecast survival profile for a single patient with competing risks
 *
 * This function forecasts the future survival time and exit cause for a patient,
 * given their observed data and conditional survival probabilities. It handles
 * right-censored and interval-censored cases.
 *
 * @param log_cond_prob_surv Array of row vectors containing log conditional 
 *        probabilities of survival for each cause and time interval.
 *        Dimensions: [n_causes, n_intervals]
 * @param exit_cause Observed exit cause (1 to n_causes)
 * @param event_time Observed event time or time of censoring
 * @param right_censored Indicator if the patient is right-censored (1 if censored, 0 otherwise)
 * @param interval_censored Number of intervals over which the observation is interval censored (0 if not interval censored)
 *
 * @return A tuple containing:
 *         1. int: The forecasted last interval before exit (survival time)
 *         2. int: Indicator if the forecast is right censored (1 if censored, 0 otherwise)
 *         3. int: The forecasted cause of exit (1 to n_causes)
 *
 * The function handles three scenarios:
 * 1. For right-censored patients, it forecasts the future survival profile
 * 2. For interval-censored patients, it draws a survival time within the censored interval
 * 3. For patients with observed exits, it returns the observed data
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

