/**
 * Data Preparation for Confirmed Response and Tumor Analysis
 *
 * This block of code prepares and validates data for a survival analysis model
 * that incorporates confirmed response times and tumor characteristics.
 */

// Calculate the maximum week for confirmed response
int max_confresp_week = max(max(confirmed_response_week), extend_max_confresp_week);

print("max_confresp_week = ", max_confresp_week);

// Arrays and vector for time interval distances between baseline hazards
array[max_confresp_week] real confresp_range;      // Time in years (assuming 12 intervals per year)
array[max_confresp_week] int confresp_range_int;   // Time in integer intervals
vector[max_confresp_week] confresp_range_vec;      // Time in intervals (zero-indexed)

// Initialize time ranges
for (i in 1:max_confresp_week) {
  confresp_range[i] = i / 12.0;
  confresp_range_int[i] = i;
  confresp_range_vec[i] = i - 1;  
}

// Prepare tumor covariate matrices
matrix[n_patients, n_tumor_covar] tumor_sum_covar;  
matrix[n_patients, n_tumor_covar] uncentered_tumor_sum_covar; // Just scaled
vector[2] tumor_sum_covar_mean;
vector[2] tumor_sum_covar_sd;

// Calculate tumor sum covariates
(tumor_sum_covar, uncentered_tumor_sum_covar, tumor_sum_covar_mean, tumor_sum_covar_sd) =  
  prepare_early_tumor_sums_covar(tumor_size, n_patient_tumors, n_measures, t_measure, n_screening_t, n_tumor_covar);  

// Calculate total confirmed response weeks and time periods
int<lower = 0> n_total_confresp_week = sum(confirmed_response_week);
int<lower = 0> n_crcr_time_periods = max_confresp_week * n_patients;

// Arrays for last unclassified response week and confirmed response cause
array[n_patients] int<lower = 0> last_unclassified_response_week;
array[n_patients] int<lower = 1, upper = n_causes> confirmed_response_cause;

// Calculate last unclassified response week and confirmed response cause for each patient
for (i in 1:n_patients) {
  last_unclassified_response_week[i] =  
    confirmed_response_week[i] - (1 - confirmed_response_censored[i]) * (confirmed_response_interval_censored[i] * (1 - crcr_ignore_interval_censoring) + 1);
  confirmed_response_cause[i] = confirmed_response[i] + 1;
}

// Create array of patient confirmed response interval positions
array[n_patients + 1] int<lower = 1> patient_conf_resp_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_confresp_week + 1);

// Calculate indices for observed and missing confirmed response values
int n_missing_confirmed_response = sum(confirmed_response_censored);  
int n_obs_confirmed_response = n_patients - n_missing_confirmed_response;  
array[n_obs_confirmed_response] int<lower = 1, upper = n_patients> obs_confirmed_response;
array[n_missing_confirmed_response] int<lower = 1, upper = n_patients> missing_confirmed_response;

(obs_confirmed_response, missing_confirmed_response) = get_mask_idx(confirmed_response_censored);

// Arrays for trial-specific objective response rate (ORR) population
array[n_trials] int<lower = 0> n_trial_orr_pop;
array[sum(orr_pop)] int<lower = 1> trial_orr_pop;

// Calculate trial-specific ORR population
{
  int trial_orr_pop_pos = 1;
  
  for (s in 1:n_trials) {
    int patient_pos = trial_patient_pos[s];
    int patient_end = trial_patient_pos[s + 1] - 1;
  
    n_trial_orr_pop[s] = sum(orr_pop[patient_pos:patient_end]);
  
    print("trial ", s, ": n = ", n_trial_patients[s], ", n_orr_pop = ", n_trial_orr_pop[s]);
  
    for (i in 1:n_trial_patients[s]) {
      if (orr_pop[patient_pos + i - 1]) {
        trial_orr_pop[trial_orr_pop_pos] = i;
        trial_orr_pop_pos += 1;
      }
    }
  }
}