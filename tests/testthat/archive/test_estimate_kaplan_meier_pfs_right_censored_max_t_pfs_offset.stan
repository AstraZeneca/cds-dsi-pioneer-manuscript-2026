functions {
  // Copy the specific estimate_kaplan_meier function we want to test
  tuple(vector, array[] int, array[] int, array[] int) estimate_kaplan_meier(array[] int pfs, array[] int right_censored, int max_t, int pfs_offset) {
    int n_pfs = size(pfs); // How many patients
    array[n_pfs] int sorted_pfs_idx = sort_indices_asc(pfs);
    int pfs_pos = 1;
    int n = n_pfs; // How many patients still haven't seen disease progression. 
    
    vector[max_t + 1] s = rep_vector(1.0, max_t + 1);
    array[max_t + 1] int at_risk = rep_array(n, max_t + 1);
    array[max_t + 1] int n_right_censored = rep_array(0, max_t + 1);
    array[max_t + 1] int n_exited = rep_array(0, max_t + 1); 
    
    for (t in 1:(max_t + 1)) {
      real prev_s = t > 1 ? s[t - 1] : 1.0;
      
      while ((n > 0) && (pfs_pos <= n_pfs) && (pfs[sorted_pfs_idx[pfs_pos]] <= t)) {
        if (t <= max_t) {
          // Remember that we define "pfs" as the last interval survived not the interval of exit.
          n_exited[t + pfs_offset] += !right_censored[sorted_pfs_idx[pfs_pos]]; 
        }
        
        n_right_censored[t] += right_censored[sorted_pfs_idx[pfs_pos]];
       
        pfs_pos += 1; 
        
      }
    
      s[t] = n > 0 ? prev_s * (n - n_exited[t]) / n : prev_s;
      at_risk[t] = n; 
      n -= n_exited[t] + n_right_censored[t];
    }
    
    return (s, at_risk, n_right_censored, n_exited); 
  }
}

data {
  int<lower=0> n_patients;
  array[n_patients] int<lower=0> pfs;
  array[n_patients] int<lower=0,upper=1> right_censored;
  int<lower=0> max_t;
  int<lower=0> pfs_offset;
}

generated quantities {
  vector[max_t + 1] km_survival;
  array[max_t + 1] int at_risk;
  array[max_t + 1] int n_right_censored;
  array[max_t + 1] int n_exited;
  
  (km_survival, at_risk, n_right_censored, n_exited) = estimate_kaplan_meier(pfs, right_censored, max_t, pfs_offset);
}

