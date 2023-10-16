functions {
  #include "util.stan"
}

data {
  int<lower = 0, upper = 1> fit_data;
  int<lower = 0, upper = 1> gen_tumor_sizes;
  
  #include "base_data.stan"
  
  matrix<lower = 0>[gen_tumor_sizes ? 0 : sum(n_patient_tumors), n_measures] tumor_size;
}

transformed data {
  array[n_measures] real all_measure_idx;
  real delta = 1e-9;
  real<lower = 0> tumor_sd = rows(tumor_size) > 1 ? sd(tumor_size) : 0;
  matrix<lower = 0>[rows(tumor_size) > 1 ? sum(n_patient_tumors) : 0, n_measures] scaled_tumor_size = tumor_size / tumor_sd; 
  
  for (m in 1:n_measures) {
    all_measure_idx[m] = m;
  }
}

parameters {
  real pop_tumor_gp_intercept;
  real<lower = 0> pop_tumor_gp_alpha;
  real<lower = 0> pop_tumor_gp_rho;
  real<lower = 0> pop_tumor_sigma;

  matrix[n_measures, sum(n_patient_tumors)] eta;  
}

transformed parameters {
  matrix[sum(n_patient_tumors), n_measures] tumor_pred; 
  
  {
    int tumor_pos = 1;
    
    for (i in 1:n_patients) {
      for (j in 1:n_patient_tumors[i]) {
        tumor_pred[tumor_pos] = calc_gp_pred(all_measure_idx, pop_tumor_gp_intercept, pop_tumor_gp_alpha, pop_tumor_gp_rho, delta, eta[, tumor_pos])'; 
        
        tumor_pos += 1;
      }
    }
  }
}

model {
  pop_tumor_gp_intercept ~ normal(0, 0.25);
  pop_tumor_gp_alpha ~ normal(0, 0.25);
  pop_tumor_gp_rho ~ inv_gamma(5, 5);
  pop_tumor_sigma ~ normal(0, 0.5);
  
  to_vector(eta) ~ std_normal();

  if (fit_data) {
    to_vector(scaled_tumor_size) ~ lognormal(to_vector(tumor_pred), pop_tumor_sigma);
  }  
}

generated quantities {
  matrix<lower = 0>[gen_tumor_sizes ? sum(n_patient_tumors) : 0, n_measures] rep_tumor_size;
  
  if (gen_tumor_sizes) {
    rep_tumor_size = to_matrix(lognormal_rng(to_vector(tumor_pred), pop_tumor_sigma), sum(n_patient_tumors), n_measures);
  }
}

