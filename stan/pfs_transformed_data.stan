/**
 * Data Preparation and Validation for Survival Analysis
 *
 * This block of code prepares and validates data for a survival analysis model,
 * specifically for a piecewise constant hazard model with censoring.
 */

// Total sum of progression-free survival (PFS) times across all patients
int<lower = 0, upper = max(pfs) * n_patients> n_total_pfs = sum(pfs);

// Arrays and vector for time interval distances between baseline hazards
array[max_all_t] real pfs_range;      // Time in years (assuming 12 intervals per year)
array[max_all_t] int pfs_range_int;   // Time in integer intervals
vector[max_all_t] pfs_range_vec;      // Time in intervals (zero-indexed)

// Initialize time ranges
for (i in 1:max_all_t) {
  pfs_range[i] = i / 12.0;
  pfs_range_int[i] = i;
  pfs_range_vec[i] = i - 1;  
}

// Array to indicate non-right-censored patients
array[n_patients] int<lower = 0, upper = 1> right_uncensored = rep_array(0, n_patients);

{
  // Validate censoring and create right_uncensored array
  for (i in 1:n_patients) {
      right_uncensored[i] = 1 - right_censored[i];
  
      // Check for invalid censoring (both interval and right censored)
      if (interval_censored[i] > 0 && right_censored[i]) {
        reject("Cannot be both interval and right censored.");
      }
  }
  
  // Print censoring summary
  print("Number of interval censored observations: ", sum(interval_censored));
  print("Number of right censored observations: ", sum(right_censored));
}

// Calculate the number of time periods for the model
// If generating PFS, use all possible time intervals; otherwise, use observed PFS plus censored observations
int<lower = 0> n_time_periods = gen_pfs ? n_patients * max_all_t : n_total_pfs + sum(right_uncensored) + sum(interval_censored);