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
  // which
  array[N_CASES, MAX_LEN] int which_out;
  // all_eq
  array[N_CASES] int all_eq_out;
  // any_eq
  array[N_CASES] int any_eq_out;
  // sort_asc
  array[N_CASES, MAX_LEN] int sort_asc_out;
  // reverse
  array[N_CASES, MAX_LEN] int reverse_out;
  // which_min
  array[N_CASES] int which_min_out;
  // which_max
  array[N_CASES] int which_max_out;
  // is_sorted
  array[N_CASES] int is_sorted_out;
  // is_unique
  array[N_CASES] int is_unique_out;
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
    // which (find indices where x == 1)
    which_out[case] = rep_array(-9999, MAX_LEN);
    if (n_x[case] > 0) {
      array[MAX_LEN] int which_tmp = which(x[case,1:n_x[case]], 1);
      for (i in 1:MAX_LEN) which_out[case,i] = which_tmp[i];
    }
    // all_eq (are all x == 1?)
    all_eq_out[case] = all_eq(x[case,1:n_x[case]], 1);
    // any_eq (is any x == 1?)
    any_eq_out[case] = any_eq(x[case,1:n_x[case]], 1);
    // sort_asc
    sort_asc_out[case] = rep_array(-9999, MAX_LEN);
    if (n_x[case] > 0) {
      array[n_x[case]] int sort_tmp = sort_asc(x[case,1:n_x[case]]);
      for (i in 1:n_x[case]) sort_asc_out[case,i] = sort_tmp[i];
    }
    // reverse
    reverse_out[case] = rep_array(-9999, MAX_LEN);
    if (n_x[case] > 0) {
      array[n_x[case]] int rev_tmp = reverse(x[case,1:n_x[case]]);
      for (i in 1:n_x[case]) reverse_out[case,i] = rev_tmp[i];
    }
    // which_min
    which_min_out[case] = which_min(x[case,1:n_x[case]]);
    // which_max
    which_max_out[case] = which_max(x[case,1:n_x[case]]);
    // is_sorted
    is_sorted_out[case] = is_sorted(x[case,1:n_x[case]]);
    // is_unique
    is_unique_out[case] = is_unique(x[case,1:n_x[case]]);
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
