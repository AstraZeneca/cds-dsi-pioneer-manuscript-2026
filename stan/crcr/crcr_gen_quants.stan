array[n_patients] matrix[max_confresp_week, n_causes] cif; // cumulative incidence function
matrix[n_patients, n_causes] prob_cause; 

(cif, prob_cause) = calc_cif(n_patients, log_crcr_cond_prob_surv, max_confresp_week);