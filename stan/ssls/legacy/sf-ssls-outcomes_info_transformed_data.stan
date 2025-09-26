array[n_pfs_timepoints] int<lower = 0> sorted_pfs_timepoints = sort_asc(pfs_timepoints);

array[n_cond_group + 1] int cond_group_pos = create_pos(cond_group_size);
