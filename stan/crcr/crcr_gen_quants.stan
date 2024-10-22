array[n_patients] matrix[max_confresp_week, n_causes] cif; // cumulative incidence function
matrix[n_patients, n_causes] prob_cause; 

(cif, prob_cause) = calc_cif(n_patients, log_crcr_cond_prob_surv, max_confresp_week);

array[n_patients] int<lower = 1> rep_confirmed_response_week;
array[n_patients] int<lower = 0, upper = 1> rep_confirmed_response_censored;
array[n_patients] int<lower = 0, upper = 2> rep_confirmed_response; // 2 happens when censored
array[n_patients] int<lower = 0, upper = 1> rep_confirmed_response_forced; // Censored entries are randomly assigned a response based off CIF 

for (i in 1:n_patients) {
  int patient_prob_pos = 1 + (i - 1) * max_confresp_week; 
  int patient_prob_end = patient_prob_pos + max_confresp_week - 1;
  
  (rep_confirmed_response_week[i], rep_confirmed_response_censored[i], rep_confirmed_response[i]) = 
    competing_risks_survival_time_rng(log_crcr_cond_prob_surv[patient_prob_pos:patient_prob_end]);
    
  rep_confirmed_response[i] -= 1;
  rep_confirmed_response_week[i] += 1 - rep_confirmed_response_censored[i];
  
  rep_confirmed_response_forced[i] = rep_confirmed_response_censored[i] ? bernoulli_rng(prob_cause[2])[1] : rep_confirmed_response[i];
}


real<lower = 0, upper = 1> rep_orr = mean(rep_confirmed_response_forced);