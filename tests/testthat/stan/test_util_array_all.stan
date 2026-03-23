functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
}
data {
  int<lower=1> N;
  int<lower=1> N_REPEAT;
  array[N] int arr_int;
  array[N] int mask;
  array[N_REPEAT] int to_repeat;
  int repeats;
  int mon;
  int first_calendar;
  int calendar_date;
  array[N] int id_arr;
  int<lower=1> N_ALL;
  array[N_ALL] int search_arr;
  int search_what;
  int n_succ;
}
generated quantities {
  // count_positive
  int cnt_pos = count_positive(arr_int);

  // which (normal: indices where mask > 0)
  array[count_positive(mask)] int which_idx = which(mask);

  // which inverse (indices where mask <= 0)
  array[N - count_positive(mask)] int which_inv_idx = which(mask, 1);

  // rep_each
  array[N_REPEAT * repeats] int rep_each_out = rep_each(to_repeat, repeats);

  // months_to_weeks
  real weeks = months_to_weeks(mon);

  // calendar_date_to_study_date (scalar overload)
  int study_date = calendar_date_to_study_date(first_calendar, calendar_date);

  // calendar_date_to_study_date (vectorised: array first_dates, scalar calendar_date)
  array[N] int study_dates_vec = calendar_date_to_study_date(arr_int, calendar_date);

  // id2idx
  array[N] int idx_arr = id2idx(id_arr);
}
