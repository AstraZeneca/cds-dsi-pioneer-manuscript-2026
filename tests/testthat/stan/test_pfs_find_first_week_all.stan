functions {
  #include "multistate.stanfunctions"
  #include "pfs.stanfunctions"
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}

// Harness tests two functions:
//   find_first_week(arr, values, min_run_length, curr_visits, forecast_time, max_all_t)
//   find_first_forecast_week(arr, values, min_run_length, forecast_time, max_all_t)
//
// N_CASES / N_FC_CASES use padded 2D arrays; case-specific lengths are passed via n_* arrays.
// Slicing to 1:0 produces a valid empty array when n_obs[case] = 0.

data {
  // ---- find_first_week cases ----
  int<lower=1> N_CASES;
  int<lower=1> MAX_ARR;
  int<lower=1> MAX_VALS;
  int<lower=1> MAX_OBS;
  int<lower=1> MAX_FORE;
  array[N_CASES] int n_arr;
  array[N_CASES] int n_vals;
  array[N_CASES] int n_obs;
  array[N_CASES] int n_fore;
  array[N_CASES] int min_run_length;
  array[N_CASES] int max_all_t;
  array[N_CASES, MAX_ARR]  int arr;
  array[N_CASES, MAX_VALS] int vals;
  array[N_CASES, MAX_OBS]  int curr_visits;
  array[N_CASES, MAX_FORE] int forecast_time;

  // ---- find_first_forecast_week cases ----
  int<lower=1> N_FC_CASES;
  int<lower=1> MAX_ARR_FC;
  int<lower=1> MAX_VALS_FC;
  int<lower=1> MAX_FORE_FC;
  array[N_FC_CASES] int n_arr_fc;
  array[N_FC_CASES] int n_vals_fc;
  array[N_FC_CASES] int n_fore_fc;
  array[N_FC_CASES] int min_run_fc;
  array[N_FC_CASES] int max_all_t_fc;
  array[N_FC_CASES, MAX_ARR_FC]  int arr_fc;
  array[N_FC_CASES, MAX_VALS_FC] int vals_fc;
  array[N_FC_CASES, MAX_FORE_FC] int forecast_time_fc;
}

generated quantities {
  // find_first_week outputs
  array[N_CASES] int out_week;
  array[N_CASES] int out_rc;

  for (case in 1:N_CASES) {
    (out_week[case], out_rc[case]) = find_first_week(
      arr[case, 1:n_arr[case]],
      vals[case, 1:n_vals[case]],
      min_run_length[case],
      curr_visits[case, 1:n_obs[case]],
      forecast_time[case, 1:n_fore[case]],
      max_all_t[case]
    );
  }

  // find_first_forecast_week outputs
  array[N_FC_CASES] int out_week_fc;
  array[N_FC_CASES] int out_rc_fc;

  for (case in 1:N_FC_CASES) {
    (out_week_fc[case], out_rc_fc[case]) = find_first_forecast_week(
      arr_fc[case, 1:n_arr_fc[case]],
      vals_fc[case, 1:n_vals_fc[case]],
      min_run_fc[case],
      forecast_time_fc[case, 1:n_fore_fc[case]],
      max_all_t_fc[case]
    );
  }
}
