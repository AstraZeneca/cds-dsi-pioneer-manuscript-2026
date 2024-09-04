vector<lower = 0>[n_bootstrap_param] bootstrap_median_pfs;
array[n_bootstrap_param] int<lower = 0, upper = n_patients> n_bootstrap_sample;
vector<lower = 0, upper = 1>[n_bootstrap_param] bootstrap_maturity_rate;

{ 
  int pfs_interval_pos = 1;
  array[n_causes] matrix[max_all_t, n_patients] mat_log_cond_prob_surv;
  
  for (i in 1:n_patients) {
    int n_intervals = max_all_t;
    int pfs_interval_end = pfs_interval_pos + n_intervals - 1;
    
    for (k in 1:n_causes) {
      mat_log_cond_prob_surv[k, , i] = log_cond_prob_surv[pfs_interval_pos:pfs_interval_end, k];
    }
    
    pfs_interval_pos = pfs_interval_end + 1;
  }
  
  for (p in 1:n_bootstrap_param) {
    bootstrap_maturity_rate[p] = 0;

    array[n_patients] int experiment_start = sort_asc(neg_binomial_2_rng(rep_vector(recruit_lambda[1], n_patients), rep_vector(recruit_phi[1], n_patients)));
    array[n_patients] int current_bootstrap_pfs;

    int i = 0;

    while ((i + 1 <= n_patients) && (experiment_start[i + 1] <= prediction_week[p])) {
      i += 1;

      int current_patient = discrete_range_rng(1, n_patients);
      int bootstrap_confirmed_response;

      if (experiment_start[i] + confirmed_response_week[current_patient] - 1 <= prediction_week[p]) {
        bootstrap_confirmed_response = confirmed_response[current_patient];
      } else {
        // Not observed yet, so let's estimate it.
        bootstrap_confirmed_response = bernoulli_rng(conf_resp_prob[current_patient, 2]);
      }

      bootstrap_maturity_rate[p] += bootstrap_confirmed_response;
      current_bootstrap_pfs[i] = survival_time_rng(mat_log_cond_prob_surv[bootstrap_confirmed_response + 1, , current_patient]).1;
    }

    n_bootstrap_sample[p] = i;

    if (n_bootstrap_sample[p] > 0) {
      bootstrap_median_pfs[p] = survival_median(current_bootstrap_pfs[:n_bootstrap_sample[p]], max_all_t).1;
      bootstrap_maturity_rate[p] /= n_bootstrap_sample[p];
    } else {
      bootstrap_median_pfs[p] = 0;
      bootstrap_maturity_rate[p] = 0;
    }
  }
}