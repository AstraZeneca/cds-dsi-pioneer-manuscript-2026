// _other_events_transformed_data.inc
// Recreated from tumor/other_events_transformed_data.stan for modular ssls build.
// Provides non-target PFS event timing and right-censor indicators restricted
// to the full patient set (legacy logic preserved).
array[n_patients] int<lower = 0> non_target_pfs; 
array[n_patients] int<lower = 0, upper = 1> non_target_right_censored = zeros_int_array(n_patients);

for (i in 1:n_patients) {
  non_target_pfs[i] = pfs[i] + interval_censored[i];
  non_target_right_censored[i] = pfs[i] >= target_pfs[i];
}
