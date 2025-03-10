int gen_pfs = 1;

#include "../base_transformed_data.stan" 
#include "../pfs_transformed_data.stan"
#include "../crcr/crcr_transformed_data.stan"

int crcr_grain_size = 83; // For reduce_sum()

array[n_patients + 1] int<lower = 1> patient_pfs_interval_pos = linspaced_int_array(n_patients + 1, 1, n_patients * max_all_t + 1);