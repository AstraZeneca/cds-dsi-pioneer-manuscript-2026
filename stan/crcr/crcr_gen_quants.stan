array[n_patients] matrix<lower = 0>[max_confresp_week, n_causes] cif; // cumulative incidence function. Couldn't add an upper constraint because in prior predict we get 1 + epsilons.
matrix<lower = 0, upper = 1>[n_patients, n_causes] prob_cause; // Approx probability of exiting to each of the competing risks

(cif, prob_cause) = calc_cif(n_patients, log_crcr_cond_prob_surv, max_confresp_week);

// Posterior predicted outcomes
array[n_patients] int<lower = 1> rep_confirmed_response_week;
array[n_patients] int<lower = 0, upper = 1> rep_confirmed_response_censored;
array[n_patients] int<lower = 0, upper = 2> rep_confirmed_response; // 2 happens when censored
array[n_patients] int<lower = 0, upper = 1> rep_confirmed_response_forced; // Censored entries are randomly assigned a response based off CIF 

// Forecast outcomes
array[n_patients] int<lower = 1> forecast_confirmed_response_week;
array[n_patients] int<lower = 0, upper = 1> forecast_confirmed_response_censored;
array[n_patients] int<lower = 0, upper = 2> forecast_confirmed_response; // 2 happens when censored
array[n_patients] int<lower = 0, upper = 1> forecast_confirmed_response_forced; // Censored entries are randomly assigned a response based off CIF 

for (i in 1:n_patients) {
  int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
  int patient_prob_end = patient_prob_pos + max_confresp_week - 1;
  
  matrix[max_confresp_week, n_causes] patient_log_crcr_cond_prob_surv = log_crcr_cond_prob_surv[patient_prob_pos:patient_prob_end];
  
  (rep_confirmed_response_week[i], rep_confirmed_response_censored[i], rep_confirmed_response[i]) = competing_risks_survival_time_rng(patient_log_crcr_cond_prob_surv);
    
  rep_confirmed_response[i] -= 1;
  rep_confirmed_response_week[i] += 1 - rep_confirmed_response_censored[i];
  rep_confirmed_response_forced[i] = rep_confirmed_response_censored[i] ? bernoulli_rng(prob_cause[2])[1] : rep_confirmed_response[i];
  
  if (confirmed_response_censored[i] || confirmed_response_interval_censored[i] > 0) {
    (forecast_confirmed_response_week[i], forecast_confirmed_response_censored[i], forecast_confirmed_response[i]) = competing_risks_survival_time_rng(
      patient_log_crcr_cond_prob_surv, confirmed_response[i] + 1, last_unclassified_response_week[i], confirmed_response_censored[i], confirmed_response_interval_censored[i]
    );
    
    forecast_confirmed_response[i] -= 1;
    forecast_confirmed_response_week[i] += 1 - forecast_confirmed_response_censored[i];
    
    forecast_confirmed_response_forced[i] = forecast_confirmed_response_censored[i] ? bernoulli_rng(prob_cause[2])[1] : forecast_confirmed_response[i];
  } else { // Use what is observed
    forecast_confirmed_response_censored[i] = 0;
    forecast_confirmed_response[i] = confirmed_response[i];
    forecast_confirmed_response_forced[i] = confirmed_response[i];
    forecast_confirmed_response_week[i] = confirmed_response_week[i];
  }
}

// Objective response rate
real<lower = 0, upper = 1> rep_orr = mean(rep_confirmed_response_forced); 
vector<lower = 0, upper = 1>[n_trials] rep_trial_orr;

real<lower = 0, upper = 1> forecast_orr = mean(forecast_confirmed_response_forced);
vector<lower = 0, upper = 1>[n_trials] forecast_trial_orr;

for (s in 1:n_trials) {
  int patient_pos = trial_patient_pos[s];
  int patient_end = trial_patient_pos[s + 1] - 1;
  
  rep_trial_orr[s] = mean(rep_confirmed_response_forced[patient_pos:patient_end]);
  forecast_trial_orr[s] = mean(forecast_confirmed_response_forced[patient_pos:patient_end]);
}