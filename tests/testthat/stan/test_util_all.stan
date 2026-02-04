functions {
  #include "util.stan"
  #include "pos.stan"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_LEN;
  // Explicit test case dimensions from R
  array[N_CASES] int<lower=1> n_patients_case;
  array[N_CASES] int<lower=1> n_tumors_case;
  array[N_CASES] int<lower=1> n_meas_sum_case;
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

  // summarize_matrix_eigenvalues
  matrix[2,2] test_mat;
  test_mat[1,1] = 2; test_mat[1,2] = 0; test_mat[2,1] = 0; test_mat[2,2] = 8;
  tuple(real, real, real) eig = summarize_matrix_eigenvalues(test_mat);
  real min_eig = eig.1;
  real max_eig = eig.2;
  real cond_num = eig.3;

  // calculate_n_missing_measures
  array[2] int n_measures_test = {2, 2};
  array[4] int t_measure_test = {1, 2, 3, 4};
  array[2] int n_patient_tumors_test = {1, 1};
  array[2] int n_missing_measures = calculate_n_missing_measures(n_measures_test, t_measure_test, n_patient_tumors_test);

  // get_mask_idx
  array[5] int mask = {0, 1, 0, 1, 1};
  tuple(array[2] int, array[3] int) mask_idx = get_mask_idx(mask);
  array[2] int idx0 = mask_idx.1;
  array[3] int idx1 = mask_idx.2;

  // unique_by_pos
  array[6] int x_test = {1, 2, 2, 3, 4, 4};
  array[3] int pos_test = {1, 4, 7};
  tuple(array[4] int, array[3] int) uniq_by_pos = unique_by_pos(x_test, pos_test);
  array[4] int uniq_vals = uniq_by_pos.1;
  array[3] int uniq_pos = uniq_by_pos.2;

  // get_idx_dict
  array[4] int idx = {2, 4, 4, 7};
  array[7] int idx_dict = get_idx_dict(idx);

  // standardize_tumor_sizes
  vector[3] tumor_size = [1.0, 2.0, 3.0]';
  tuple(real, real, vector[3]) std_tumor = standardize_tumor_sizes(tumor_size);
  real mean_tumor = std_tumor.1;
  real sd_tumor = std_tumor.2;
  vector[3] std_vals = std_tumor.3;

  for (case in 1:N_CASES) {

    // Use explicit test case dimensions from R
    int n_patients = n_patients_case[case];
    int n_tumors = n_tumors_case[case];
    int n_meas_sum = n_meas_sum_case[case];
    array[n_tumors] int n_meas;
    int t_meas_pos = 1;
    array[n_meas_sum] int t_meas;
    for (i in 1:n_tumors) {
      n_meas[i] = n_measures[case, i];
      for (j in 1:n_meas[i]) {
        t_meas[t_meas_pos] = t_measure[case, t_meas_pos];
        t_meas_pos += 1;
      }
    }
    array[n_patients] int n_pat_tum;
    for (i in 1:n_patients) n_pat_tum[i] = n_patient_tumors[case, i];
    for (i in 1:MAX_LEN) max_t_out[case, i] = 0;
    if (n_patients > 0 && n_meas_sum > 0) {
      array[n_patients] int max_t = get_max_t(t_meas, n_meas, n_pat_tum);

      for (i in 1:n_patients) max_t_out[case, i] = max_t[i];
    }
    // num_unique
    int n_unique = num_unique(x[case,1:n_x[case]]);
    num_unique_out[case] = n_unique;

    // unique
    n_unique_out[case] = n_unique;
    for (i in 1:MAX_LEN) unique_out[case,i] = 0;
    if (n_x[case] > 0 && n_unique > 0) {
      array[n_unique] int uniq = unique(x[case,1:n_x[case]]);

      for (i in 1:n_unique) unique_out[case,i] = uniq[i];
    }
    // find_first
    find_first_out[case] = find_first(all[case,1:n_all[case]], what[case,1:n_what[case]], n_succ[case]);

  }
}
