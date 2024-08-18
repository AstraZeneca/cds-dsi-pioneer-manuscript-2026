data {
  int<lower = 0, upper = 1> fit_data;
  
  int<lower = 0> n_trial_sim; 
  
  // This is the data that is shared with the tumor model 
  #include "base_data.stan"
   
  array[n_patients] int<lower = 1> experiment_start_week; // Week 1 is the first week of the experiment
}

parameters {
  vector<lower = 0>[n_trials] alpha; 
  vector<lower = 0>[n_trials] beta; 
  vector<lower = 0>[n_trials] lambda;
  vector<lower = 0>[n_trials] phi;
}

model {
  // Priors
  
  alpha ~ normal(50, 5); // Prior for the shape parameter
  beta ~ normal(1, 1); // Prior for the rate parameter
  lambda ~ gamma(alpha, beta); 
  phi ~ normal(0, 50);
  
  // Likelihood
  
  if (fit_data) {
    // experiment_start_week ~ poisson(lambda[patient_trial]);  
    // experiment_start_week ~ neg_binomial(alpha[patient_trial], beta[patient_trial]);  
    experiment_start_week ~ neg_binomial_2(lambda[patient_trial], phi[patient_trial]);  
  }
}

generated quantities {
  vector<lower = 0>[n_trials] recruit_sd = sqrt(lambda + (lambda^2) ./ phi);
  array[n_trials, n_trial_sim] int<lower = 0> sim_experiment_start_week;
 
  if (n_trial_sim > 0) {
    for (s in 1:n_trials) {
      // sim_experiment_start_week[s] = poisson_rng(rep_vector(lambda[s], n_trial_sim));
      // sim_experiment_start_week[s] = neg_binomial_rng(rep_vector(alpha[s], n_trial_sim), rep_vector(beta[s], n_trial_sim));
      sim_experiment_start_week[s] = neg_binomial_2_rng(rep_vector(lambda[s], n_trial_sim), rep_vector(phi[s], n_trial_sim));
    }
  } 
}