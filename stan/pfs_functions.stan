
/**
 * Calculate the marginal probability of disease progression at every interval.
 *
 * This function computes the log marginal probability of disease progression
 * (or "exit") for each time interval, given the log conditional probabilities
 * of survival.
 *
 * @param log_cond_prob_surv A row vector containing the log conditional 
 *                           probabilities of survival for each time interval.
 * @return A row vector of log marginal probabilities of disease progression 
 *         for each time interval.
 */
row_vector calculate_log_marginal_exit_prob(row_vector log_cond_prob_surv) {
  int T = num_elements(log_cond_prob_surv); 
  row_vector[T] log_marginal_exit_prob;  // Output vector

  for (t in 1:T) {
    // Calculate the log probability of exit at time t
    log_marginal_exit_prob[t] = log1m_exp(log_cond_prob_surv[t]);
  
    if (t > 1) {
      // For t > 1, add the sum of log conditional survival probabilities
      // up to t-1 to account for survival up to the previous interval
      log_marginal_exit_prob[t] += sum(log_cond_prob_surv[1:(t - 1)]);
    }
  }
  
  return(log_marginal_exit_prob);
}

/**
 * Calculate the log-likelihood for a piecewise constant hazard model.
 *
 * This function computes the log-likelihood for each patient in a survival analysis
 * using a piecewise constant hazard model. It handles various censoring scenarios
 * including right censoring and interval censoring.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param exit_event Array indicating the type of exit event for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Array of matrices containing log conditional survival probabilities for each exit type and patient
 * @param start_from Array of starting intervals for each patient
 * @param end_at Array of ending intervals for each patient
 *
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring,  
  array[] matrix log_cond_prob_surv, array[] int start_from, array[] int end_at
) {
  int n_exit_types = size(log_cond_prob_surv), n_patients = rows(log_cond_prob_surv[1]);
  vector[n_patients] lp = zeros_vector(n_patients);
  
  // Iterate over all patients
  for (i in 1:n_patients) {
    int interval_pos = max(0, start_from[i]);
    int interval_end = min(end_at[i], last_surv_week[i]);
    
    // Given the value of end_at, this patient might become right censored.
    int effective_right_censored = right_censored[i] || (end_at[i] < last_surv_week[i] + (1 - ignore_interval_censoring) * interval_censored[i] + 1);
    int curr_interval_censored = ignore_interval_censoring || effective_right_censored ? 0 : interval_censored[i];
  
    // Check for invalid censoring scenario
    if (right_censored[i] && interval_censored[i] > 0) {
      fatal_error("Interval censoring not allowed with right censored observations. Patient ", i, ".");
    }
  
    // Calculate log-likelihood for known survival part
    if (interval_pos <= interval_end) {
      for (k in 1:n_exit_types) {
        lp[i] += sum(log_cond_prob_surv[k, i, interval_pos:interval_end]);
      }
    }
  
    // Handle interval censoring
    if ((interval_pos <= interval_end) || (interval_end + curr_interval_censored + 1 >= interval_pos)) {
      int old_interval_end = interval_end, old_interval_censored = curr_interval_censored;
      interval_end = max(old_interval_end, interval_pos - 1);
      curr_interval_censored = max(0, curr_interval_censored - (interval_end - old_interval_end));  
  
      // Create a mixture of all possible true intervals of exit for interval censoring
      vector[curr_interval_censored + 1] ic_mix_lp = zeros_vector(curr_interval_censored + 1);  
  
      for (c in 0:curr_interval_censored) {
        if (c > 0) {  
          for (k in 1:n_exit_types) {
            ic_mix_lp[c + 1] += sum(log_cond_prob_surv[k, i, (interval_end + 1):(interval_end + c)]);
          }
        }
  
        if (!effective_right_censored) {  
          ic_mix_lp[c + 1] += log1m_exp(log_cond_prob_surv[exit_event[i], i, interval_end + c + 1]);
        }
      }
  
      // Add log-likelihood for interval censoring
      lp[i] += curr_interval_censored > 0 ? log_sum_exp(ic_mix_lp) - log(curr_interval_censored + 1) : ic_mix_lp[1];  
    }
  }
  
  return lp;
}

/**
 * Overloaded versions of calc_pch_loglik function for different parameter combinations.
 * These functions provide convenience wrappers around the main calc_pch_loglik function,
 * allowing it to be called with fewer parameters by using default values or derived parameters.
 */

/**
 * Calculate log-likelihood with specified start times but default end times.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param exit_event Array indicating the type of exit event for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Array of matrices containing log conditional survival probabilities for each exit type and patient
 * @param start_from Array of starting intervals for each patient
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring,  
  array[] matrix log_cond_prob_surv, array[] int start_from
) {
  int max_all_t = cols(log_cond_prob_surv[1]);
  
  return calc_pch_loglik(
    last_surv_week, exit_event, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv, start_from, rep_array(max_all_t, size(start_from))
  );
}  

/**
 * Calculate log-likelihood with default start times.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param exit_event Array indicating the type of exit event for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Array of matrices containing log conditional survival probabilities for each exit type and patient
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, array[] matrix log_cond_prob_surv
) {
  return calc_pch_loglik(last_surv_week, exit_event, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_surv, ones_int_array(size(last_surv_week)));
}

/**
 * Calculate log-likelihood for a single exit type with default start times.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Matrix containing log conditional survival probabilities for each patient
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, matrix log_cond_prob_surv
) {
  int n_patients = size(last_surv_week);
  
  return calc_pch_loglik(last_surv_week, ones_int_array(n_patients), right_censored, interval_censored, ignore_interval_censoring, { log_cond_prob_surv }, ones_int_array(n_patients));
}

/**
 * Calculate log-likelihood for a single exit type with specified start times.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Matrix containing log conditional survival probabilities for each patient
 * @param start_from Array of starting intervals for each patient
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, matrix log_cond_prob_surv, array[] int start_from
) {
  int n_patients = size(last_surv_week);
  
  return calc_pch_loglik(last_surv_week, ones_int_array(n_patients), right_censored, interval_censored, ignore_interval_censoring, { log_cond_prob_surv }, start_from);
}

/**
 * Calculate log-likelihood for a single exit type with specified start and end times.
 *
 * @param last_surv_week Array of last survival weeks for each patient
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0)
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1
 * @param log_cond_prob_surv Matrix containing log conditional survival probabilities for each patient
 * @param start_from Array of starting intervals for each patient
 * @param end_at Array of ending intervals for each patient
 * @return A vector of log-likelihoods, one for each patient
 */
vector calc_pch_loglik(
  array[] int last_surv_week, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, matrix log_cond_prob_surv,  
  array[] int start_from, array[] int end_at
) {
  int n_patients = size(last_surv_week);
  
  return calc_pch_loglik(last_surv_week, ones_int_array(n_patients), right_censored, interval_censored, ignore_interval_censoring, { log_cond_prob_surv }, start_from, end_at);
}

/**
 * Calculate the piecewise-constant proportional hazard log-likelihood if the log conditional survival probabilities are stored in a single vector.
 *
 * @param pfs Array of observed number of weeks without disease progression for each patient.
 * @param right_censored Array indicating if each patient is right-censored (1) or not (0).
 * @param interval_censored Array indicating the number of weeks over which we have interval censoring for each patient.
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1.
 * @param log_cond_prob_surv Vector of log conditional probabilities of survival at every interval.
 * @param max_all_t The latest week assessment is done across all the data.
 * @param start_from Array of weeks in which the first post-treatment assessment was done for each patient.
 *
 * @return Vector of patient-level log-likelihoods.
 *
 * Note: This function differs from the pch_lpmf() function in that it returns the vector of 
 * patient-level log-likelihoods rather than a single log-likelihood value for all patients.
 */
vector calc_pch_loglik(
  array[] int pfs,  
  array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring,  
  vector log_cond_prob_surv, int max_all_t, array[] int start_from  
) {
  int n_patients = size(pfs);
  
  return calc_pch_loglik(
    pfs, ones_int_array(n_patients), right_censored, interval_censored, ignore_interval_censoring, { to_matrix(log_cond_prob_surv, max_all_t, n_patients) }, start_from
  );
}

/**
 * Piecewise Constant Hazard Log Probability Mass Function
 *
 * These functions provide an easy-to-use Stan distribution for piecewise constant hazard models.
 * They calculate the total log-likelihood by summing individual log-probabilities.
 */

/**
 * @param y Array of observed number of weeks without disease progression.
 * @param right_censored Array indicating if each observation is right-censored (1) or not (0).
 * @param interval_censored Array of number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1.
 * @param log_cond_prob_progress Vector of log conditional probabilities of disease progression at every interval.
 * @param max_all_t The latest week assessment is done across all the data.
 * @param start_from Array of weeks in which the first post-treatment assessment was done.
 * @return Total log-likelihood summed across all observations.
 */
real pch_lpmf(
  array[] int y,  
  array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, vector log_cond_prob_progress, int max_all_t, array[] int start_from
) {
  return sum(calc_pch_loglik(y, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, max_all_t, start_from));
}

/**
 * @param y Array of observed number of weeks without disease progression.
 * @param right_censored Array indicating if each observation is right-censored (1) or not (0).
 * @param interval_censored Array of number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1.
 * @param log_cond_prob_progress Array of matrices containing log conditional probabilities of disease progression.
 * @param start_from Array of weeks in which the first post-treatment assessment was done.
 * @return Total log-likelihood summed across all observations.
 */
real pch_lpmf(
  array[] int y,  
  array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, array[] matrix log_cond_prob_progress, array[] int start_from
) {
  return sum(calc_pch_loglik(y, ones_int_array(size(y)), right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, start_from));
}

/**
 * @param y Array of observed number of weeks without disease progression.
 * @param exit_event Array indicating the type of exit event for each observation.
 * @param right_censored Array indicating if each observation is right-censored (1) or not (0).
 * @param interval_censored Array of number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1.
 * @param log_cond_prob_progress Array of matrices containing log conditional probabilities of disease progression.
 * @param start_from Array of weeks in which the first post-treatment assessment was done.
 * @param end_at Array of weeks in which the last assessment was done.
 * @return Total log-likelihood summed across all observations.
 */
real pch_lpmf(
  array[] int y,  
  array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, array[] matrix log_cond_prob_progress,  
  array[] int start_from, array[] int end_at
) {
  return sum(calc_pch_loglik(y, exit_event, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, start_from, end_at));
}

/**
 * @param y Array of observed number of weeks without disease progression.
 * @param exit_event Array indicating the type of exit event for each observation.
 * @param right_censored Array indicating if each observation is right-censored (1) or not (0).
 * @param interval_censored Array of number of weeks over which we have interval censoring.
 * @param ignore_interval_censoring Integer flag to ignore interval censoring if set to 1.
 * @param log_cond_prob_progress Array of matrices containing log conditional probabilities of disease progression.
 * @param start_from Array of weeks in which the first post-treatment assessment was done.
 * @return Total log-likelihood summed across all observations.
 */
real pch_lpmf(
  array[] int y,  
  array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, array[] matrix log_cond_prob_progress, array[] int start_from
) {
  return sum(calc_pch_loglik(y, exit_event, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, start_from));
}

real pch_lpmf(
  array[] int y, 
  array[] int exit_event, array[] int right_censored, array[] int interval_censored, int ignore_interval_censoring, array[] matrix log_cond_prob_progress
) {
  return sum(calc_pch_loglik(y, exit_event, right_censored, interval_censored, ignore_interval_censoring, log_cond_prob_progress, ones_int_array(size(y))));
}
/**
 * Calculate the quantiles of the given sample PFS (Progression-Free Survival).
 *
 * This function computes specified quantiles of survival times, accounting for right censoring.
 * It's particularly useful for calculating median PFS and other percentiles of interest.
 *
 * @param surv_time Array of patient survival times.
 * @param last_surv_time Interval beyond which we assume right censoring. This may not affect lower quantiles.
 * @param p Vector of percentiles to calculate (e.g., [0.25, 0.5, 0.75] for quartiles).
 * @return A tuple containing:
 *         1. A vector of calculated quantiles.
 *         2. An array indicating which quantiles couldn't be calculated due to right censoring (1 if not calculable, 0 if calculable).
 */
tuple(vector, array[] int) survival_quantiles(array[] int surv_time, int last_surv_time, vector p) {
  int N = size(surv_time);
  int P = num_elements(p); // Number of percentiles
  vector[P] quantiles;
  array[P] int cannot_calculate;
  array[N] int sorted_surv_time = sort_asc(surv_time);
  
  // Calculate each quantile
  for (j in 1:P) {
      int k = 1;
  
      // We can't use the floor() function because it returns real and to_int() expected "data real".
      while (p[j] >= k * (1.0 / N)) {
        k += 1;
      }
  
      if (sorted_surv_time[max(k - 1, 1)] > last_surv_time) {
        quantiles[j] = 0.0;
        cannot_calculate[j] = 1;
      } else {
        real pos = p[j] * (N - 1) + 1;
        real d = pos - (k - 1);
  
        quantiles[j] = sorted_surv_time[max(k - 1, 1)] + d * (sorted_surv_time[k] - sorted_surv_time[max(k - 1, 1)]);
        cannot_calculate[j] = 0;
      }
  }
  
  return (quantiles, cannot_calculate);
}

/**
 * Calculate median PFS (Progression-Free Survival)
 *
 * This function is a specialized version of survival_quantiles that calculates only the median (50th percentile).
 * It accounts for right censoring in the data.
 *
 * @param surv_time Array of patient survival times.
 * @param last_surv_time Interval beyond which we assume right censoring. This may not affect the median calculation in many cases.
 * @return A tuple containing:
 *         1. The calculated median survival time.
 *         2. An integer flag indicating whether the median could be calculated (0 if calculable, 1 if not calculable due to right censoring).
 */
tuple(real, int) survival_median(array[] int surv_time, int last_surv_time) {
  vector[1] q;
  array[1] int c;
  (q, c) = survival_quantiles(surv_time, last_surv_time, [ 0.5 ]');
  
  return(q[1], c[1]);
}

/** Calculate the proportion of patients who survived beyond time time n (PFSn).
 * 
 * @param surv_time Patient survival times
 * @param n
 * @return Proportion surviving >= n
 */
real calc_pfs_n(array[] int surv_time, real n) {
  int n_patients = size(surv_time);
  array[n_patients] int sorted_surv_time = sort_desc(surv_time);
  int pfs_n = 0;
  
  while (pfs_n < n_patients && sorted_surv_time[pfs_n + 1] >= n) {
    pfs_n += 1;
  }
  
  return 1.0 * pfs_n / n_patients;
}

/**
 * Generate survival times given log conditional probabilities of survival
 *
 * This function simulates a survival time based on the provided log conditional 
 * probabilities of survival for each time interval.
 *
 * @param log_cond_prob_surv Row vector of log conditional probabilities of survival at each interval
 * @return A tuple containing:
 *         1. The simulated survival time (integer)
 *         2. A censoring indicator (1 if right-censored, 0 otherwise)
 */
tuple(int, int) survival_time_rng(row_vector log_cond_prob_surv) {
  int n_intervals = num_elements(log_cond_prob_surv);
  
  int survival_time = 0;
  
  row_vector[n_intervals] log_odds_surv = log_cond_prob_surv - log1m_exp(log_cond_prob_surv);
  
  while (survival_time < n_intervals && bernoulli_logit_rng(log_odds_surv[survival_time + 1])) {
    survival_time += 1;
  }
  
  int censored = survival_time >= n_intervals;  
  
  return(survival_time, censored);
}

/**
 * Forecast survival time based on observed data and log conditional probabilities
 *
 * This function forecasts the survival time for a patient given their observed 
 * survival time, censoring status, and the log conditional probabilities of survival.
 *
 * @param log_cond_prob_surv Row vector of log conditional probabilities of survival
 * @param obs_surv_time Observed survival time
 * @param right_censored Indicator for right censoring (1 if right-censored, 0 otherwise)
 * @param interval_censored Number of intervals over which the observation is interval censored (0 if not interval censored)
 * @return A tuple containing:
 *         1. The forecasted survival time (integer)
 *         2. An indicator for whether the forecast is right-censored (1 if right-censored, 0 otherwise)
 */
tuple(int, int) survival_time_rng(row_vector log_cond_prob_surv, int obs_surv_time, int right_censored, int interval_censored) {
  int n = num_elements(log_cond_prob_surv);
  int survival_time = obs_surv_time, forecast_right_censored = right_censored;
  
  if (right_censored) {
    // Check that we are not at the highest interval; we don't have any estimated probability after that.  
    if (obs_surv_time < n) {
      (survival_time, forecast_right_censored) = survival_time_rng(log_cond_prob_surv[(obs_surv_time + 1):]);
      survival_time += obs_surv_time;
    }
  } else if (interval_censored > 0) {
    survival_time += interval_censored_survival_time_rng(log_cond_prob_surv[(obs_surv_time + 1):(obs_surv_time + interval_censored + 1)]);
  }
  
  return(survival_time, forecast_right_censored);
}


tuple(int, int) survival_time_rng(row_vector log_cond_prob_surv, int obs_surv_time) {
  return survival_time_rng(log_cond_prob_surv, obs_surv_time, 1, 0);
}

/**
 * Generate a survival time within the range of interval censored intervals
 *
 * Have to be careful using this; it is not in the log space. I don't want to touch it right now because I don't want to break something critical.
 *
 * This function simulates a survival time for an interval censored observation,
 * given the log conditional probabilities of survival within the censored interval.
 *
 * @param ic_log_cond_prob_surv Row vector of log conditional probabilities of survival 
 *                              at each interval within the interval censored range
 * @return An integer representing the time after the left bound of the interval censored range
 */  
int interval_censored_survival_time_rng(row_vector ic_log_cond_prob_surv) {
  int n = num_elements(ic_log_cond_prob_surv);
  row_vector[n] marginal_prob_exit = exp(cumulative_sum(append_col(0, ic_log_cond_prob_surv))[:n] + log1m_exp(ic_log_cond_prob_surv));

  marginal_prob_exit /= sum(marginal_prob_exit);

  return categorical_rng(marginal_prob_exit') - 1;
}

/**
 * Estimate Kaplan-Meier survival function aggregated over all patients
 *
 * This function calculates the Kaplan-Meier estimate of the survival function S(t) = Pr[T > t],
 * where t ranges from 0 to max_t. It accounts for right-censoring and multiple causes of exit.
 *
 * @param last_surv Array of integers representing the last observed week that was progression-free for each patient
 * @param cause Array of integers indicating the cause of exit for each patient
 * @param right_censored Array of integers (0 or 1) indicating whether each patient was right-censored
 * @param max_t Integer specifying the last interval to report Kaplan-Meier results
 *
 * @return A tuple containing four elements:
 *         1. vector: Proportion surviving at each time point (S(t))
 *         2. array[] int: Number of patients at risk at each time point
 *         3. array[] int: Number of right-censored patients at each time point
 *         4. array[,] int: Number of patients who exited due to each cause at each time point
 *
 * Note: The returned arrays and vector are of length max_t + 1, with index 1 corresponding to t=0.
 */
tuple(vector, array[] int, array[] int, array[,] int) estimate_kaplan_meier(array[] int last_surv, array[] int cause, array[] int right_censored, int max_t) {
  int n_patients = size(last_surv); // How many patients
  int n_causes = size(cause);
  array[n_patients] int sorted_last_surv_idx = sort_indices_asc(last_surv);
  int last_surv_pos = 1;
  int n = n_patients;

  vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
  array[max_t + 1] int at_risk = rep_array(n, max_t + 1);
  array[max_t + 1] int n_right_censored = rep_array(0, max_t + 1);
  array[max_t + 1, n_causes] int n_exited;

  // For each time interval in 0..max_t see how many patients exited and calculate proportion surviving.

  for (t in 0:max_t) {
    n_exited[t + 1] = rep_array(0, n_causes);
    real prev_s = t > 0 ? s[t] : 1.0;

    while (
      (n > 0) &&  
      (last_surv_pos <= n_patients) &&  
      (right_censored[sorted_last_surv_idx[last_surv_pos]] || (last_surv[sorted_last_surv_idx[last_surv_pos]] <= t))
    ) {
      n_exited[t + 1, cause[sorted_last_surv_idx[last_surv_pos]]] += !right_censored[sorted_last_surv_idx[last_surv_pos]];
      n_right_censored[t + 1] += right_censored[sorted_last_surv_idx[last_surv_pos]];

      last_surv_pos += 1;
    }

    s[t + 1] = n > 0 ? prev_s * (n - sum(n_exited[t + 1])) / n : prev_s;
    at_risk[t + 1] = n;
    n -= sum(n_exited[t + 1]) + n_right_censored[t + 1];
  }

  return (s, at_risk, n_right_censored, n_exited);
}

/** Survival aggregated over all patients, S(t) = Pr[T > t], t \in {0,..., N} 
 * For each time interval in 0..max_t see how many patiented exited and calculate proportion surviving.
 * 
 * This should probably be a wrapper for the above function. Won't touch it now.
 *
 * @param pfs The last observed week that was progression-free
 * @param right_censored Right censoring per patient
 * @param max_t The last interval to report Kaplan-Meier results
 * @return (Proportion surviving, Number at risk, Number right censored, Number for whom disease progressed) for each week
 */
tuple(vector, array[] int, array[] int, array[] int) estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t) {
  int n_pfs = size(pfs); // How many patients
  array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
  int pfs_pos = 1;
  int n = n_pfs; // How many patients still haven't seen disease progression. 
  
  vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
  array[max_t + 1] int at_risk = rep_array(n, max_t + 1);
  array[max_t + 1] int n_right_censored = rep_array(0, max_t + 1);
  array[max_t + 1] int n_exited = rep_array(0, max_t + 1); 
  
  for (t in 1:(max_t + 1)) {
    real prev_s = t > 1 ? s[t - 1] : 1.0;
    
    while ((n > 0) && (pfs_pos <= n_pfs) && (pfs[sorted_pfs_idx[pfs_pos]] <= t)) {
      if (t <= max_t) {
        // Remember that we define "pfs" as the last interval survived not the interval of exit.
        n_exited[t + 1] += !right_censored[sorted_pfs_idx[pfs_pos]]; 
      }
      
      n_right_censored[t] += right_censored[sorted_pfs_idx[pfs_pos]];
     
      pfs_pos += 1; 
      
    }
  
    s[t] = n > 0 ? prev_s * (n - n_exited[t]) / n : prev_s;
    at_risk[t] = n; 
    n -= n_exited[t] + n_right_censored[t];
  }
  
  return (s, at_risk, n_right_censored, n_exited); 
}  

/**
 * Calculate the Concordance Index (C-index) for survival data
 *
 * The C-index measures the discriminative power of a risk score in predicting survival times.
 * It represents the probability that, for a pair of randomly chosen patients, the patient with 
 * the higher risk score has a shorter survival time.
 *
 * @param pfs Array of progression-free survival times for each patient
 * @param right_censored Array indicating whether each patient is right-censored (1) or not (0)
 * @param risk_score Vector of risk scores for each patient
 * @return The calculated C-index as a real number between 0 and 1
 */
real calc_c_index(array[] int pfs, array[] int right_censored, vector risk_score) {
  int n_patients = size(pfs);
  int n_concord = 0;
  int n_ranked = 0;
  
  for (i in 1:n_patients) {
    for (j in 1:n_patients) {
      if (i != j && pfs[i] > pfs[j] && !right_censored[j]) {
        n_ranked += 1;
        n_concord += risk_score[i] < risk_score[j];  
      }
    }
  }
  
  return 1.0 * n_concord / n_ranked;
}

/**
 * Calculate the C-index for competing risks survival data
 *
 * This version of the C-index calculation handles competing risks and incorporates 
 * information about confirmed responses.
 *
 * @param pfs Array of progression-free survival times for each patient
 * @param right_censored Array indicating whether each patient is right-censored (1) or not (0)
 * @param confirmed_response Array indicating the type of confirmed response for each patient
 * @param confirmed_response_censored Array indicating whether the confirmed response is censored for each patient
 * @param log_risk_score Array of vectors containing log risk scores for each cause and patient
 * @param log_last_cif Array of vectors containing log cumulative incidence function values for each cause and patient
 * @return The calculated C-index as a real number between 0 and 1
 */
real calc_c_index(
  array[] int pfs, array[] int right_censored,  
  array[] int confirmed_response, array[] int confirmed_response_censored,  
  array[] vector log_risk_score, array[] vector log_last_cif  
) {
  int n_patients = size(pfs);
  vector[n_patients] risk_score;
  
  for (i in 1:n_patients) {
    risk_score[i] =  
      confirmed_response_censored[i] ?  
      log_sum_exp(log_last_cif[1, i] + log_risk_score[1, i], log_last_cif[2, i] + log_risk_score[2, i]) - log_sum_exp(log_last_cif[1, i], log_last_cif[2, i]) :
      log_risk_score[confirmed_response[i] + 1, i];
  }
  
  return calc_c_index(pfs, right_censored, risk_score);
}

/**
 * Calculate the Administrative Brier Score for survival data with competing risks
 *
 * The Brier Score measures the accuracy of probabilistic predictions in survival analysis.
 * This function calculates the Brier Score at multiple time points, accounting for 
 * administrative censoring and competing risks.
 *
 * @param pfs Array of progression-free survival times for each patient
 * @param admin_right_censored_week Array of administrative right censoring times for each patient
 * @param interval_censored Array indicating the length of interval censoring for each patient
 * @param cause Array indicating the cause of event for each patient
 * @param cause_right_censored Array indicating whether the cause is right-censored for each patient
 * @param ignore_interval_censoring Flag to ignore interval censoring if set to 1
 * @param log_last_cif Array of vectors containing log cumulative incidence function values for each cause and patient
 * @param log_cond_prob_surv Array of matrices containing log conditional survival probabilities for each cause and patient
 * @return A matrix of Brier Scores, where rows represent patients and columns represent time points
 */
matrix calc_admin_brier_score(
  array[] int pfs, array[] int admin_right_censored_week, array[] int interval_censored,  
  array[] int cause, array[] int cause_right_censored,  
  int ignore_interval_censoring,  
  array[] vector log_last_cif, array[] matrix log_cond_prob_surv
) {
  int n_patients = rows(log_cond_prob_surv[1]);
  int T = cols(log_cond_prob_surv[1]);
  int n_causes = size(log_cond_prob_surv);
 
  if (n_causes > 2) {
    fatal_error("Only supports a max of two causes.");
  }
  
  matrix[n_patients, T] brier_score_t = rep_matrix(0, n_patients, T);
  
  for (t in 1:T) {
    real brier_scale = 0;
    
    for (i in 1:n_patients) {
      if (admin_right_censored_week[i] >= t) {
        brier_scale += 1;
        
        int curr_interval_censored = ignore_interval_censoring ? 0 : interval_censored[i]; 
        real observed_event = 1.0 * min(max(0, pfs[i] + curr_interval_censored + 1 - t), curr_interval_censored + 1) / (curr_interval_censored + 1);
        
        real log_risk_score; // Prob[T > t]
       
        if (n_causes > 1) {
          if (cause_right_censored[i]) {
            log_risk_score = log_sum_exp(
              log_last_cif[1, i] + log_sum_exp(calculate_log_marginal_exit_prob(log_cond_prob_surv[1, i, :t])), 
              log_last_cif[2, i] + log_sum_exp(calculate_log_marginal_exit_prob(log_cond_prob_surv[2, i, :t]))) - 
              log_sum_exp(log_last_cif[1, i], log_last_cif[2, i]);
          } else {
            log_risk_score = log_sum_exp(calculate_log_marginal_exit_prob(log_cond_prob_surv[cause[i], i, :t]));
          }
        } else {
          log_risk_score = log_sum_exp(calculate_log_marginal_exit_prob(log_cond_prob_surv[1, i, :t]));
        }
        
        brier_score_t[i, t] = square((1 - exp(log_risk_score)) - observed_event);
      }
    }
   
    if (brier_scale > 0) {
      brier_score_t[, t] /= brier_scale; 
    } else { // No more patients with admin censoring after t: set the loss to zero.
      break;
    }
    
  }
  
  return brier_score_t;
}

/** Convert a ragged vector tumor size measures to a matrix aligned by measurement week
 *
 * @param tumor_size Tumor sizes
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_screening_t Number of observed pre-screening assessments per tumor
 * @param max_measures Number of tumor sizes to include in output, irrespective of when they are actually observed.
 * @return (Aligned matrix, Actual t used for the columns)
 */
tuple(matrix, array[,] int) create_aligned_matrix(vector tumor_size, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures) {
  int n_tumors = size(n_measures);
  int n_total_measures = sum(n_measures);
  array[n_total_measures] int t_sort_ind = sort_indices_asc(t_measure);
  array[n_total_measures] int within_tumor_ind;
  
  int last_within_ind = 0;
  int last_t = min(t_measure) - 1; 
  
  for(n in 1:n_total_measures) {
    if (t_measure[t_sort_ind[n]] != last_t) {
      last_within_ind += 1;
    }
    
    last_t = t_measure[t_sort_ind[n]];
    within_tumor_ind[t_sort_ind[n]] = last_within_ind;
  }

  int n_screening_col = max(n_screening_t);
  int n_aligned_col = max(max_measures, last_within_ind); 
  matrix[n_tumors, n_aligned_col] aligned_matrix = rep_matrix(0, n_tumors, n_aligned_col);
  array[n_tumors, n_aligned_col] int aligned_t = rep_array(0, n_tumors, n_aligned_col);
  
  int tumor_pos = 1;
  
  for (j in 1:n_tumors) {
    int tumor_end = tumor_pos + n_measures[j] - 1;
    
    aligned_matrix[j, within_tumor_ind[tumor_pos:tumor_end]] = tumor_size[tumor_pos:tumor_end]'; 
    aligned_t[j, within_tumor_ind[tumor_pos:tumor_end]] = t_measure[tumor_pos:tumor_end]; 
    
    tumor_pos = tumor_end + 1;
  }
  
  return (aligned_matrix[, n_screening_col:(max_measures + n_screening_col - 1)], aligned_t[, n_screening_col:(max_measures + n_screening_col - 1)]);
}

/** Create (n_patients * n_tumors) x max_measures matrix of each tumor's covariates from assessment **order** (not t) = 1, 2, ...., max_measures 
 * What this function actually returns is the last pre-screening measure and the (max_measures - 1) succeeding measures. The second output in the tuple
 * provides the intervals (weeks) of these max_measures columns.
 * 
 * @param tumor_size Tumor sizes
 * @param n_patient_tumors Array with the number of tumors per patient.
 * @param n_measures The number of assessments per tumor.
 * @param t_measure The week each assessment was done.
 * @param n_screening_t Number of observed pre-screening assessments per tumor
 * @param max_measures Number of tumor sizes to include in output, irrespective of when they are actually observed.
 * @return (Design matrix, Actual t used for the columns)
 */
tuple(matrix, array[,] int) prepare_early_tumors_design_matrix(
  vector tumor_size, array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
  matrix[sum(n_patient_tumors), max_measures] tumor_covar = rep_matrix(0, sum(n_patient_tumors), max_measures);
  array[sum(n_patient_tumors), max_measures] int tumor_covar_t = rep_array(min(t_measure) - 1, sum(n_patient_tumors), max_measures);
  int n_patients = size(n_patient_tumors);
  int tumor_pos = 1;
  int tumor_size_pos = 1;
  int covar_pos = 1;
  
  for (i in 1:n_patients) {
    int n_current_tumors = n_patient_tumors[i];
    int covar_end = covar_pos + n_current_tumors - 1; 
    int tumor_end = tumor_pos + n_current_tumors - 1;
    int n_patient_measures = sum(n_measures[tumor_pos:tumor_end]);
    int tumor_size_end = tumor_size_pos + n_patient_measures - 1;
   
    (tumor_covar[covar_pos:covar_end], tumor_covar_t[covar_pos:covar_end]) = create_aligned_matrix(
      tumor_size[tumor_size_pos:tumor_size_end], n_measures[tumor_pos:tumor_end], t_measure[tumor_size_pos:tumor_size_end], n_screening_t, max_measures
    );
    
    covar_pos = covar_end + 1;
    tumor_pos = tumor_end + 1;
    tumor_size_pos = tumor_size_end + 1;
  }
  
  return (tumor_covar, tumor_covar_t);
}

tuple(matrix, matrix, vector, vector) prepare_early_tumor_sums_covar(
  vector tumor_size, 
  array[] int n_patient_tumors, array[] int n_measures, array[] int t_measure, array[] int n_screening_t, int max_measures
) {
  int n_tumors = sum(n_patient_tumors);
  int n_patients = size(n_patient_tumors);
  matrix[n_tumors, 2] tumor_covar; 
  array[n_tumors, 2] int tumor_covar_t; 
  matrix[n_patients, 2] tumor_sum_covar; 
  matrix[n_patients, 2] uncentered_tumor_sum_covar; 
  
  (tumor_covar, tumor_covar_t) = prepare_early_tumors_design_matrix(tumor_size, n_patient_tumors, n_measures, t_measure, n_screening_t, max_measures);
  
  {
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      int tumor_end = tumor_pos + n_patient_tumors[i] - 1;
      
      tumor_sum_covar[i] = ones_row_vector(n_patient_tumors[i]) * tumor_covar[tumor_pos:tumor_end];
      
      tumor_pos = tumor_end + 1;
    }
  }
  
  vector[2] tumor_sum_covar_mean;
  vector[2] tumor_sum_covar_sd;
  
  for (c in 1:2) {
    uncentered_tumor_sum_covar[, c] = tumor_sum_covar[, c];
    
    (tumor_sum_covar_mean[c], tumor_sum_covar_sd[c], tumor_sum_covar[, c]) = standardize_tumor_sizes(tumor_sum_covar[, c]);
    
    uncentered_tumor_sum_covar[, c] /= tumor_sum_covar_sd[c];
  }
  
  return (tumor_sum_covar, uncentered_tumor_sum_covar, tumor_sum_covar_mean, tumor_sum_covar_sd);
}