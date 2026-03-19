functions {
  #include "pos.stanfunctions"
  #include "util.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  int<lower=1> n_cases;
  int<lower=1> max_len;
  array[n_cases] int n_trajectory;
  array[n_cases] int n_screening;
  array[n_cases, max_len] real sld_trajectory;
  array[n_cases] real pre_nadir;
}



generated quantities {
  array[n_cases, max_len] int result;
  array[n_cases, max_len] int result_with_nadir;
  for (i in 1:n_cases) {
    int n_treat = n_trajectory[i] - n_screening[i];
    vector[n_trajectory[i]] sld = to_vector(sld_trajectory[i, 1:n_trajectory[i]]);
    array[max_len] int tmp = rep_array(-1, max_len);
    array[max_len] int tmp_nadir = rep_array(-1, max_len);
    if (n_treat > 0) {
      array[n_treat] int res = calculate_target_recist(sld, n_screening[i]);
      array[n_treat] int res_nadir = calculate_target_recist(sld, pre_nadir[i], n_screening[i]);
      for (j in 1:n_treat) {
        tmp[j] = res[j];
        tmp_nadir[j] = res_nadir[j];
      }
    }
    for (j in 1:max_len) {
      result[i, j] = tmp[j];
      result_with_nadir[i, j] = tmp_nadir[j];
    }
  }
}

