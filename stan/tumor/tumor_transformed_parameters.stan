// Given priors on actual tumor mean and SD calculate the corresponding lognormal ones
real<lower = 0> pop_tumor_gp_alpha = sqrt(log((tumor_sd / tumor_mean)^2 + 1));
real pop_tumor_gp_intercept = log(tumor_mean) - pop_tumor_gp_alpha^2 / 2.0;

vector[use_tumor_model ? sum(n_patient_tumors) : 0] tumor_gp_intercept; 

if (use_tumor_model) {
  tumor_gp_intercept = rep_vector(pop_tumor_gp_intercept, sum(n_patient_tumors));
}

if (use_tumor_model && multilevel_patient) {
  int tumor_pos = 1;
  
  for (i in 1:n_patients) {
    for (j in 1:n_patient_tumors[i]) {
      tumor_gp_intercept[tumor_pos] += patient_tumor_gp_intercept_effect[i];
      
      tumor_pos += 1;
    }
  }
}

if (use_tumor_model && multilevel_tumor) {
  tumor_gp_intercept += tumor_gp_intercept_effect;
}