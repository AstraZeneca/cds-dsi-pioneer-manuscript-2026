array[n_tumor_separate_trials] vector[model_all_measures ? max_t_width : n_pop_unique_visits] latent_pop_tumor_gp;

for (s in 1:n_tumor_separate_trials) {
  latent_pop_tumor_gp[s] = calc_gp_pred(
    model_all_measures ? all_measure_t : all_measure_t[pop_unique_visits_idx], 0, pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], delta, latent_pop_tumor_gp_eta[s]
  ); 
}

// Given priors on actual tumor mean and SD calculate the corresponding lognormal ones
// vector<lower = 0>[n_tumor_separate_trials] pop_tumor_sigma = sqrt(log((tumor_sd ./ tumor_mean)^2 + 1));
// vector[n_tumor_separate_trials] pop_tumor_gp_intercept = log(tumor_mean) - pop_tumor_sigma^2 / 2.0;

// vector[n_patients] tumor_gp_intercept;

// array[n_trials] matrix[max(patient_max_t_width), max(patient_max_t_width)] trial_tumor_gp_cov = rep_array(diag_matrix(rep_vector(1, max(patient_max_t_width))), n_trials);
// array[n_trials] cov_matrix[max(patient_max_t_width)] trial_tumor_gp_cov;

// for (s in 1:n_trials) {
//   int curr_patient_pos = trial_patient_pos[s];
//   int curr_patient_end = trial_patient_pos[s + 1] - 1;
//   
//   // tumor_gp_intercept[curr_patient_pos:curr_patient_end] = rep_vector(pop_tumor_gp_intercept[min(s, n_tumor_separate_trials)], n_trial_patients[s]);
//   
//   // if (add_trial_level_tumor_gp) {
//   //   pop_tumor_gp_intercept[curr_patient_pos:curr_patient_end] += trial_tumor_gp_intercept_effect[s];
//   // }
//   // 
//   // if (add_patient_level_tumor_gp) {
//   //   pop_tumor_gp_intercept[curr_patient_pos:curr_patient_end] += patient_tumor_gp_intercept_effect[s];
//   // }
//   
//   // Population Ks
//   if (s <= n_tumor_separate_trials) {
//     trial_tumor_gp_cov[s] = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], 1e-9);
//     // trial_tumor_gp_cov[s] = calc_gp_vcov(all_measure_t, pop_tumor_gp_alpha[s], pop_tumor_gp_rho[s], pop_tumor_sigma[s]^2);
//   } else {
//     trial_tumor_gp_cov[s] = trial_tumor_gp_cov[1];
//   }
//   
//   // Multilevel stuff
//   
//   for (i in curr_patient_pos:curr_patient_end) {
//     int curr_patient_visits_pos = patient_visit_pos[i];
//     int curr_patient_visits_end = patient_visit_pos[i + 1] - 1;
//   }
// }
// 
