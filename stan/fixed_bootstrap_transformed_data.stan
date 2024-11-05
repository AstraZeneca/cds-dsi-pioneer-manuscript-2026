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