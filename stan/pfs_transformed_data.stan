int<lower = 0, upper = max(pfs) * n_patients> n_total_pfs = sum(pfs);

// These are used for the time interval distance between baseline hazards
array[max_all_t] real pfs_range;
array[max_all_t] int pfs_range_int;
vector[max_all_t] pfs_range_vec;

for (i in 1:max_all_t) {
  pfs_range[i] = i / 12.0;
  pfs_range_int[i] = i;
  pfs_range_vec[i] = i - 1; 
}

array[n_patients] int<lower = 0, upper = 1> right_uncensored = rep_array(0, n_patients);

{
  tuple(array[n_patients] int, array[n_patients] int) censoring_res = identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure);
  
  for (i in 1:n_patients) {
      right_uncensored[i] = 1 - right_censored[i];
      
      if (interval_censored[i] > 0 && right_censored[i]) {
        reject("Cannot be both interval and right censored.");
      }
  }
    
  print("Number of interval censored observations: ", sum(interval_censored));
  print("Number of right censored observations: ", sum(right_censored));
}

// If generating PFS we need to calculate probs for all possible time intervals, otherwise only up to observed PFS. 
int<lower = 0> n_time_periods = gen_pfs ? n_patients * max_all_t : n_total_pfs + sum(right_uncensored) + sum(interval_censored);
  