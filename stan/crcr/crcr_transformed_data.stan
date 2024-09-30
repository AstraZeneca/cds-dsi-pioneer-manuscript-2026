int n_causes = 2;

int max_confresp_week = max(confirmed_response_week);

print("max_confresp_week = ", max_confresp_week);

// These are used for the time interval distance between baseline hazards
array[max_confresp_week] real confresp_range;
array[max_confresp_week] int confresp_range_int;
vector[max_confresp_week] confresp_range_vec;

for (i in 1:max_confresp_week) {
  confresp_range[i] = i / 12.0;
  confresp_range_int[i] = i;
  confresp_range_vec[i] = i - 1; 
}

matrix[n_patients, n_tumor_covar] tumor_sum_covar; 
matrix[n_patients, n_tumor_covar] uncentered_tumor_sum_covar; // Just scaled
vector[2] tumor_sum_covar_mean;
vector[2] tumor_sum_covar_sd;

(tumor_sum_covar, uncentered_tumor_sum_covar, tumor_sum_covar_mean, tumor_sum_covar_sd) = 
  prepare_early_tumor_sums_covar(tumor_size, n_patient_tumors, n_measures, t_measure, n_screening_t, n_tumor_covar); 
  
int<lower = 0> n_total_confresp_week = sum(confirmed_response_week);
int<lower = 0> n_crcr_time_periods = max_confresp_week * n_patients;

// This is the number of weeks a patient "survived" without response classification
array[n_patients] int<lower = 0> last_unclassified_response_week;
array[n_patients] int<lower = 1, upper = n_causes> confirmed_response_cause;

for (i in 1:n_patients) {
  last_unclassified_response_week[i] = 
    confirmed_response_week[i] - (1 - confirmed_response_censored[i]) * (confirmed_response_interval_censored[i] * (1 - ignore_interval_censoring) + 1);
  confirmed_response_cause[i] = confirmed_response[i] + 1;
}

// Not using last_unclassified_response_week 
array[n_patients] int<lower = 0, upper = 1> early_confirmed_response_censored = confirmed_response_censored;

array[n_patients + 1] int<lower = 1> patient_conf_resp_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_confresp_week + 1);
