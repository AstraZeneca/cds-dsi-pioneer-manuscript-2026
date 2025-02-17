vector[n_patients] log_odds_confirmed_response; 

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
  array[n_causes] row_vector[max_confresp_week] patient_log_crcr_cond_prob_surv = log_crcr_cond_prob_surv[, i];
  
  (rep_confirmed_response_week[i], rep_confirmed_response_censored[i], rep_confirmed_response[i]) = competing_risks_survival_time_rng(patient_log_crcr_cond_prob_surv);
  
  log_odds_confirmed_response[i] = log_cif[2, i, max_confresp_week] - log_cif[1, i, max_confresp_week];
    
  rep_confirmed_response[i] -= 1;
  rep_confirmed_response_week[i] += 1 - rep_confirmed_response_censored[i];
  rep_confirmed_response_forced[i] = rep_confirmed_response_censored[i] ? bernoulli_logit_rng(log_odds_confirmed_response[i]) : rep_confirmed_response[i];
  
  if (confirmed_response_censored[i] || confirmed_response_interval_censored[i] > 0) {
    (forecast_confirmed_response_week[i], forecast_confirmed_response_censored[i], forecast_confirmed_response[i]) = competing_risks_survival_time_rng(
      patient_log_crcr_cond_prob_surv, confirmed_response[i] + 1, last_unclassified_response_week[i], confirmed_response_censored[i], confirmed_response_interval_censored[i]
    );
    
    forecast_confirmed_response[i] -= 1;
    forecast_confirmed_response_week[i] += 1 - forecast_confirmed_response_censored[i];
    
    forecast_confirmed_response_forced[i] = forecast_confirmed_response_censored[i] ? bernoulli_logit_rng(log_odds_confirmed_response[i]) : forecast_confirmed_response[i];
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
vector<lower = 0, upper = 1>[n_trials] forecast_trial_subpop_orr;

{
  int orr_pop_pos = 1;
  
  for (s in 1:n_trials) {
    int patient_pos = trial_patient_pos[s];
    int patient_end = trial_patient_pos[s + 1] - 1;
    
    int orr_pop_end = orr_pop_pos + n_trial_orr_pop[s] - 1; 
    
    rep_trial_orr[s] = mean(rep_confirmed_response_forced[patient_pos:patient_end]);
    forecast_trial_orr[s] = mean(forecast_confirmed_response_forced[patient_pos:patient_end]);
    forecast_trial_subpop_orr[s] = mean(forecast_confirmed_response_forced[patient_pos:patient_end][trial_orr_pop[orr_pop_pos:orr_pop_end]]);
    
    orr_pop_pos = orr_pop_end + 1;
  }
}

// Impute confirmed response status if needed 
array[n_patients] int<lower = 0, upper = 1> sim_confirmed_response = confirmed_response;
sim_confirmed_response[missing_confirmed_response] = bernoulli_logit_rng(log_odds_confirmed_response[missing_confirmed_response]); 
 