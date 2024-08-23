data {
  int<lower = 0, upper = 1> fit_data;
  
  int<lower = 0> n_trial_sim; 
  
  // This is the data that is shared with the tumor model 
  #include "../base_data.stan"
}

parameters {
  #include "recruit_parameters.stan"
}

model {
  #include "recruit_priors.stan"

  // Likelihood
  
  if (fit_data) {
    experiment_start_week ~ neg_binomial_2(recruit_lambda[patient_trial], recruit_phi[patient_trial]);  
  }
}

generated quantities {
  vector<lower = 0>[n_trials] recruit_sd = sqrt(recruit_lambda + (recruit_lambda^2) ./ recruit_phi);
  array[n_trials, n_trial_sim] int<lower = 0> sim_experiment_start_week;
 
  if (n_trial_sim > 0) {
    for (s in 1:n_trials) {
      sim_experiment_start_week[s] = neg_binomial_2_rng(rep_vector(recruit_lambda[s], n_trial_sim), rep_vector(recruit_phi[s], n_trial_sim));
    }
  } 
}
