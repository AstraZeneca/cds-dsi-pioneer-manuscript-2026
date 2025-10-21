int<lower = 0> n_pfs_timepoints;
array[n_pfs_timepoints] int<lower = 0> pfs_timepoints; // In months

// Conditioning groups or strata to estimate outcomes for a particular covar 
int<lower = 0> n_cond_group;
array[n_cond_group] int<lower = 1> cond_group_size;
array[sum(cond_group_size)] int <lower = 1, upper = n_patients> cond_group;
