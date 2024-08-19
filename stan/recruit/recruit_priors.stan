recruit_alpha ~ normal(50, 5); // Prior for the shape parameter
recruit_beta ~ normal(1, 1); // Prior for the rate parameter
recruit_lambda ~ gamma(recruit_alpha, recruit_beta); 
recruit_phi ~ normal(0, 50);
  