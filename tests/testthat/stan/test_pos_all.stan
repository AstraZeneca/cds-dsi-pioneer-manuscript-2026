functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
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
    int n = n_groups[case];
    int n_flat_case = n_flat[case];
    array[MAX_GROUPS] int n_x_full = group_sizes[case];
    array[MAX_SIZE] int x_full = flat_data[case];
    array[n] int n_x;
    array[n_flat_case] int x;
    for (i in 1:n) n_x[i] = n_x_full[i];
    for (i in 1:n_flat_case) x[i] = x_full[i];
    array[n+1] int pos = create_pos(n_x);

    // create_pos
    for (i in 1:(n+1)) pos_out[case, i] = pos[i];
    // get_pos_size
    for (i in 1:n) pos_size_out[case, i] = get_pos_size(pos, i);
    // For all group-wise outputs, only compute for groups with size > 0, else sentinel
    for (i in 1:n) {
      if (get_pos_size(pos, i) > 0) {
        max_out[case, i] = get_max_pos(x, pos, i);
        min_out[case, i] = get_min_pos(x, pos, i);
        max_idx_out[case, i] = get_max_idx(x, pos)[i];
        min_pos_out[case, i] = get_min_pos(x, pos, i);
        max_pos_out[case, i] = get_max_pos(x, pos, i);
        last_int_out[case, i] = get_int(x, pos, i, get_pos_size(pos, i));
      } else {
        max_out[case, i] = -9999;
        min_out[case, i] = -9999;
        max_idx_out[case, i] = -9999;
        min_pos_out[case, i] = -9999;
        max_pos_out[case, i] = -9999;
        last_int_out[case, i] = -9999;
      }
    }
  }
}
