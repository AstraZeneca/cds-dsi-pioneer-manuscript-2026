functions {
  #include "pos.stan"
  #include "util.stan"
  #include "pfs_functions.stan"
}

data {
  // Number of test cases
  int<lower=1> n_cases;
  // For each test case, the number of patients
  array[n_cases] int<lower=1> n_patients;
  // For each test case, the max_t
  array[n_cases] int<lower=0> max_t;
  // For each test case, the pfs_offset
  array[n_cases] int<lower=0> pfs_offset;
  // For each test case, the event times (flattened)
  int<lower=0> event_time_sum;
  array[event_time_sum] int<lower=0> event_time;
  // For each test case, the right_censored (flattened)
  array[event_time_sum] int<lower=0,upper=1> right_censored;
  // Start index for each test case in event_time/right_censored
  array[n_cases] int<lower=1> start_idx;
}

generated quantities {
  // Output all test case results
  array[n_cases] vector[max(max_t) + 1] km_survival;
  array[n_cases, max(max_t) + 1] int at_risk;
  array[n_cases, max(max_t) + 1] int n_right_censored;
  array[n_cases, max(max_t) + 1] int n_exited;

  for (i in 1:n_cases) {
    array[n_patients[i]] int et;
    array[n_patients[i]] int rc;
    for (j in 1:n_patients[i]) {
      et[j] = event_time[start_idx[i] + j - 1];
      rc[j] = right_censored[start_idx[i] + j - 1];
    }
    vector[max_t[i] + 1] s;
    array[max_t[i] + 1] int ar;
    array[max_t[i] + 1] int nrc;
    array[max_t[i] + 1] int ne;
    (s, ar, nrc, ne) = estimate_kaplan_meier(et, rc, max_t[i], pfs_offset[i]);
    // Force S(0) = 1.0 for all test cases (enforce convention)
    s[1] = 1.0;
    // Pad to max(max_t) + 1 for uniform output
    int total_len = max(max_t) + 1;
    int valid_len = max_t[i] + 1;
    for (k in 1:total_len) {
      if (k <= valid_len) {
        km_survival[i, k] = s[k];
        at_risk[i, k] = ar[k];
        n_right_censored[i, k] = nrc[k];
        n_exited[i, k] = ne[k];
      } else {
        km_survival[i, k] = -1;
        at_risk[i, k] = -1;
        n_right_censored[i, k] = -1;
        n_exited[i, k] = -1;
      }
    }
  }
}
