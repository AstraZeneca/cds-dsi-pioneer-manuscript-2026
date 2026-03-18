functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "lfo.stanfunctions"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_PATIENTS;
  int<lower=1> MAX_CUTOFFS;
  array[N_CASES] int<lower=1> n_patients;
  array[N_CASES] int<lower=1> n_cutoffs;
  array[N_CASES, MAX_PATIENTS] int sorted_last_visit_calendar_day;
  array[N_CASES, MAX_CUTOFFS]  int cutoff_calendar_day;
}

generated quantities {
  array[N_CASES, MAX_CUTOFFS] int oos_idx_out = rep_array(0, N_CASES, MAX_CUTOFFS);

  for (case in 1:N_CASES) {
    int np = n_patients[case];
    int nc = n_cutoffs[case];
    array[nc] int result = get_oos_patients_idx(
      sorted_last_visit_calendar_day[case, 1:np],
      cutoff_calendar_day[case, 1:nc]
    );
    for (c in 1:nc) {
      oos_idx_out[case, c] = result[c];
    }
  }
}
