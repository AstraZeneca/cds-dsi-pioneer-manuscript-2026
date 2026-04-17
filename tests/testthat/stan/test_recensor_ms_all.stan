functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
}

data {
  int<lower=1> N_cases;
  int<lower=1> max_n_patients;
  // Per-case patient data (rectangularized)
  array[N_cases] int<lower=1> n_patients;
  array[N_cases, max_n_patients] int ms_final_state;
  array[N_cases, max_n_patients] int ms_time_01;
  array[N_cases, max_n_patients] int ms_censored_01;
  array[N_cases, max_n_patients] int ms_time_02;
  array[N_cases, max_n_patients] int ms_time_12;
  array[N_cases, max_n_patients] int ms_time_03;
  array[N_cases, max_n_patients] int ms_time_32;
  array[N_cases, max_n_patients] int ms_os_event_12;
  array[N_cases, max_n_patients] int interval_censored;
  array[N_cases, max_n_patients] int ms_prog_deterministic;
  array[N_cases, max_n_patients] int cutoff_week;
}

generated quantities {
  // Output arrays for each case
  array[N_cases, max_n_patients] int out_final_state = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_time_01 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_censored_01 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_time_02 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_time_12 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_time_03 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_time_32 = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_ic = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_prog_det = rep_array(0, N_cases, max_n_patients);
  array[N_cases, max_n_patients] int out_ic_gap_01 = rep_array(0, N_cases, max_n_patients);

  for (c in 1:N_cases) {
    int np = n_patients[c];
    array[np] int r_fs; array[np] int r_t01; array[np] int r_c01;
    array[np] int r_t02; array[np] int r_t12; array[np] int r_t03;
    array[np] int r_t32; array[np] int r_ic; array[np] int r_pd;
    array[np] int r_ig;

    (r_fs, r_t01, r_c01, r_t02, r_t12, r_t03, r_t32, r_ic, r_pd, r_ig) =
      recensor_ms_at_cutoff(
        ms_final_state[c, 1:np], ms_time_01[c, 1:np], ms_censored_01[c, 1:np],
        ms_time_02[c, 1:np], ms_time_12[c, 1:np], ms_time_03[c, 1:np],
        ms_time_32[c, 1:np], ms_os_event_12[c, 1:np], interval_censored[c, 1:np],
        ms_prog_deterministic[c, 1:np], cutoff_week[c, 1:np]);

    for (i in 1:np) {
      out_final_state[c, i] = r_fs[i];
      out_time_01[c, i] = r_t01[i];
      out_censored_01[c, i] = r_c01[i];
      out_time_02[c, i] = r_t02[i];
      out_time_12[c, i] = r_t12[i];
      out_time_03[c, i] = r_t03[i];
      out_time_32[c, i] = r_t32[i];
      out_ic[c, i] = r_ic[i];
      out_prog_det[c, i] = r_pd[i];
      out_ic_gap_01[c, i] = r_ig[i];
    }
  }
}
