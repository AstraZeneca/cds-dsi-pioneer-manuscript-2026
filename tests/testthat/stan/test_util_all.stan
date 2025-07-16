functions {
  #include "util.stan"
  #include "pos.stan"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_LEN;
  // get_max_t
  array[N_CASES, MAX_LEN] int t_measure;
  array[N_CASES, MAX_LEN] int n_measures;
  array[N_CASES, MAX_LEN] int n_patient_tumors;
  // num_unique/unique
  array[N_CASES, MAX_LEN] int x;
  array[N_CASES] int n_x;
  // find_first
  array[N_CASES, MAX_LEN] int all;
  array[N_CASES, MAX_LEN] int what;
  array[N_CASES] int n_all;
  array[N_CASES] int n_what;
  array[N_CASES] int n_succ;
}

generated quantities {
  // Only test functions that exist in util.stan
  // get_max_t
  array[N_CASES, MAX_LEN] int max_t_out;
  // num_unique
  array[N_CASES] int num_unique_out;
  // unique
  array[N_CASES, MAX_LEN] int unique_out;
  array[N_CASES] int n_unique_out;
  // find_first
  array[N_CASES] int find_first_out;

  for (case in 1:N_CASES) {
    // get_max_t
    max_t_out[case] = rep_array(-9999, MAX_LEN);
    if (sum(n_patient_tumors[case]) > 0 && sum(n_measures[case]) > 0) {
      int n_patients = 0;
      for (i in 1:MAX_LEN) if (n_patient_tumors[case, i] > 0) n_patients += 1;
      int n_tumors = 0;
      for (i in 1:MAX_LEN) if (n_measures[case, i] > 0) n_tumors += 1;
      array[n_tumors] int n_meas;
      int n_meas_sum = 0;
      for (i in 1:n_tumors) {
        n_meas[i] = n_measures[case, i];
        n_meas_sum += n_meas[i];
      }
      array[n_meas_sum] int t_meas;
      for (i in 1:n_meas_sum) t_meas[i] = t_measure[case, i];
      array[n_patients] int n_pat_tum;
      for (i in 1:n_patients) n_pat_tum[i] = n_patient_tumors[case, i];
      array[n_patients] int max_t = get_max_t(t_meas, n_meas, n_pat_tum);
      for (i in 1:n_patients) max_t_out[case, i] = max_t[i];
    }
    // num_unique
    num_unique_out[case] = num_unique(x[case,1:n_x[case]]);
    // unique
    n_unique_out[case] = num_unique(x[case,1:n_x[case]]);
    unique_out[case] = rep_array(-9999, MAX_LEN);
    if (n_x[case] > 0) {
      array[n_unique_out[case]] int uniq = unique(x[case,1:n_x[case]]);
      for (i in 1:n_unique_out[case]) unique_out[case,i] = uniq[i];
    }
    // find_first
    find_first_out[case] = find_first(all[case,1:n_all[case]], what[case,1:n_what[case]], n_succ[case]);
  }
}
