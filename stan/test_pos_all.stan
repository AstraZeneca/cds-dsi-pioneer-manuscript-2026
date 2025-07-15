functions {
  #include "util.stan"
  #include "pos.stan"
}

data {
  int<lower=1> N_CASES;
  int<lower=1> MAX_GROUPS;
  int<lower=1> MAX_SIZE;
  array[N_CASES, MAX_GROUPS] int group_sizes;
  array[N_CASES, MAX_SIZE] int flat_data;
  array[N_CASES] int n_groups;
  array[N_CASES] int n_flat;
}

generated quantities {
  // Output arrays for each function under test
  array[N_CASES, MAX_GROUPS+1] int pos_out;
  array[N_CASES, MAX_GROUPS] int pos_size_out;
  array[N_CASES, MAX_GROUPS] int max_out;
  array[N_CASES, MAX_GROUPS] int min_out;
  array[N_CASES, MAX_GROUPS] int max_idx_out;
  array[N_CASES, MAX_GROUPS] int min_pos_out;
  array[N_CASES, MAX_GROUPS] int max_pos_out;
  array[N_CASES, MAX_GROUPS] int last_int_out;

  for (case in 1:N_CASES) {
    array[MAX_GROUPS] int n_x = group_sizes[case];
    array[MAX_GROUPS+1] int pos = create_pos(n_x);
    array[MAX_SIZE] int x = flat_data[case];
    int n = n_groups[case];
    int n_flat_case = n_flat[case];

    // create_pos
    for (i in 1:(n+1)) pos_out[case, i] = pos[i];
    // get_pos_size
    for (i in 1:n) pos_size_out[case, i] = get_pos_size(pos, i);
    // get_max, get_min, get_max_idx
    max_out[case, 1:n] = get_max(x[1:n_flat_case], pos[1:(n+1)]);
    min_out[case, 1:n] = get_min_pos(x[1:n_flat_case], pos[1:(n+1)]);
    max_idx_out[case, 1:n] = get_max_idx(x[1:n_flat_case], pos[1:(n+1)]);
    min_pos_out[case, 1:n] = get_min_pos(x[1:n_flat_case], pos[1:(n+1)]);
    max_pos_out[case, 1:n] = get_max_pos(x[1:n_flat_case], pos[1:(n+1)]);
    // get_last_int
    for (i in 1:n) {
      if (get_pos_size(pos, i) > 0) {
        last_int_out[case, i] = get_last_int(x[1:n_flat_case], pos[1:(n+1)], i);
      } else {
        last_int_out[case, i] = -9999; // sentinel for empty group
      }
    }
  }
}
