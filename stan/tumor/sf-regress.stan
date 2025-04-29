functions {
  #include "../util.stan"
  #include "../pos.stan"
}  

data {
  #include "../base_data.stan"
}

transformed data {
  int sf_rep_T = 0;
  
  #include "../base_transformed_data.stan"
  #include "tumor_transformed_data.stan"
  #include "sf-transformed_data.stan"
}

parameters {
  // Population-level parameters
  real<lower=0> pop_d;                       // Population regression rate
  real<lower=0> pop_g;                       // Population growth rate
  real<lower=0> sigma;                       // Measurement error
  
  // // Patient-specific parameters
  // array[n_patients] real<lower=0> d;         // Patient-specific regression rates
  // array[n_patients] real<lower=0> g;         // Patient-specific growth rates
  // 
  // // Hierarchical model parameters
  // real<lower=0> tau_d;                       // Between-patient SD for d
  // real<lower=0> tau_g;                       // Between-patient SD for g
}

transformed parameters {
  vector<lower = 0>[n_patients] d = rep_vector(pop_d, n_patients);         // Patient-specific regression rates
  vector<lower = 0>[n_patients] g = rep_vector(pop_g, n_patients);         // Patient-specific growth rates
}

model {
  // Priors for population parameters
  // pop_d ~ normal(0, 1);
  // pop_g ~ normal(0, 1);
  // sigma ~ normal(0, 0.1);
  pop_d ~ normal(1.5, 0.5);                 // Based on typical regression rates
  pop_g ~ normal(0.3, 0.2);                 // Based on typical growth rates
  sigma ~ normal(0.1, 0.05);                // Measurement precision 
  
  // // Priors for hierarchical SDs
  // tau_d ~ normal(0, 0.5);
  // tau_g ~ normal(0, 0.5);
  // 
  // // Patient-specific parameters come from population distribution
  // for (i in 1:n_patients) {
  //   d[i] ~ normal(pop_d, tau_d);
  //   g[i] ~ normal(pop_g, tau_g);
  // }
  
  // Likelihood
  for (i in 1:n_patients) {
    int visit_pos, visit_end;
    (visit_pos, visit_end) = get_pos(patient_visit_pos, i);
    
    for (v_idx in visit_pos:visit_end) {
      int t = t_patient_visits[v_idx];
      
      real mu = exp(-d[i] * t) + exp(g[i] * t) - 1.0;
      
      if (normalized_sld[v_idx] >= lod) {
        normalized_sld[v_idx] ~ normal(mu, sigma);
      } else {
        // target += normal_lcdf(lod | mu, sigma);  // Left-censored observation
      }
    }
  }
}

generated quantities {
  // // Generate predictions for each observation
  // array[total_obs] real y_pred;
  // 
  // for (i in 1:total_obs) {
  //   int p = patient_id[i];
  //   y_pred[i] = exp(-d[p] * t[i]) + exp(g[p] * t[i]) - 1.0;
  // }
  // 
  // // Calculate time to nadir and nadir value for each patient
  // array[n_patients] real t_min;
  // array[n_patients] real min_val;
  // 
  // for (i in 1:n_patients) {
  //   // Time to nadir from equation (5) in paper: t_min = ln(d/g)/(d - g)
  //   if (d[i] > g[i]) {  // Only calculate if mathematically valid
  //     t_min[i] = log(d[i]/g[i])/(d[i] - g[i]);
  //     min_val[i] = exp(-d[i] * t_min[i]) + exp(g[i] * t_min[i]) - 1.0;
  //   } else {
  //     t_min[i] = 0;     // No nadir if growth rate exceeds regression rate
  //     min_val[i] = 0;
  //   }
  // }
}
