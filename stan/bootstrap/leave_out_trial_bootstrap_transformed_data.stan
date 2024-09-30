int<lower = 0, upper = n_patients> n_training_patients = n_patients - (leave_out_trial > 0 ? n_trial_patients[leave_out_trial] : 0);
int<lower = 0> n_training_crcr_intervals = n_training_patients * max_confresp_week;
array[n_training_patients] int<lower = 1, upper = n_patients> training_patients; 
array[n_training_crcr_intervals] int<lower = 1, upper = n_patients * max_confresp_week> training_crcr_intervals; 

int<lower = 0, upper = n_patients> n_testing_patients = leave_out_trial > 0 ? n_trial_patients[leave_out_trial] : 0;
int<lower = 0> n_testing_crcr_intervals = n_testing_patients * max_confresp_week;
array[n_testing_patients] int<lower = 1, upper = n_patients> testing_patients; 
array[n_testing_crcr_intervals] int<lower = 1, upper = n_patients * max_confresp_week> testing_crcr_intervals; 

print("n_patients = ", n_patients, ", n_testing_patients = ", n_testing_patients);
print("n_testing_crcr_intervals = ", n_testing_crcr_intervals);

{
  int training_patient_pos = 1;
  int testing_patient_pos = 1;
  int training_interval_pos = 1;
  int testing_interval_pos = 1;
  
  for (s in 1:n_trials) {
    int n_curr_patients = n_trial_patients[s];
    int n_curr_intervals = n_curr_patients * max_confresp_week; 
    
    int first_patient = s > 1 ? sum(n_trial_patients[:(s - 1)]) + 1 : 1;
    int last_patient = first_patient + n_trial_patients[s] - 1;
    
    int first_interval = s > 1 ? sum(n_trial_patients[:(s - 1)]) * max_confresp_week + 1 : 1;
    int last_interval = first_interval + n_trial_patients[s] * max_confresp_week - 1;
    
    if (s != leave_out_trial) {
      int training_patient_end = training_patient_pos + n_curr_patients - 1;
      int training_interval_end = training_interval_pos + n_curr_intervals - 1;
      
      training_patients[training_patient_pos:training_patient_end] = linspaced_int_array(n_trial_patients[s], first_patient, last_patient); 
    
      training_crcr_intervals[training_interval_pos:training_interval_end] = 
        linspaced_int_array(n_trial_patients[s] * max_confresp_week, first_interval, last_interval); 
      
      training_patient_pos = training_patient_end + 1; 
      training_interval_pos = training_interval_end + 1;
    } else {
      int testing_patient_end = testing_patient_pos + n_curr_patients - 1;
      int testing_interval_end = testing_interval_pos + n_curr_intervals - 1;
      
      testing_patients[testing_patient_pos:testing_patient_end] = linspaced_int_array(n_trial_patients[s], first_patient, last_patient); 
      
      testing_crcr_intervals[testing_interval_pos:testing_interval_end] = 
        linspaced_int_array(n_trial_patients[s] * max_confresp_week, first_interval, last_interval); 
      
      testing_patient_pos = testing_patient_end + 1; 
      testing_interval_pos = testing_interval_end + 1;
    }
  }
}

array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_bs_sample_idx;
array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_bs_start_calendar_week; 
array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_bs_cr_mature_calendar_week; 
array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_bs_pfs_mature_calendar_week; 
array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_cr_mature_sorted_idx;
array[leave_out_trial > 0 ? n_fixed_bootstrap_samples : 0, n_bootstrap_sample] int fixed_pfs_mature_sorted_idx;

if (leave_out_trial > 0) {
  for (f in 1:n_fixed_bootstrap_samples) {
    (fixed_bs_sample_idx[f], 
     fixed_bs_start_calendar_week[f], 
     fixed_cr_mature_sorted_idx[f], 
     fixed_bs_cr_mature_calendar_week[f], 
     fixed_pfs_mature_sorted_idx[f], 
     fixed_bs_pfs_mature_calendar_week[f]) = get_bootstrap_sample_rng(
      n_bootstrap_sample, leave_out_trial, trial_patient_pos, sorted_experiment_start_week, confirmed_response_week, pfs, right_censored
    );
  }
}