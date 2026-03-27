functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "multistate.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "_hierarchy_data.stan"
  #include "_visit_data.stan"
  #include "_full_model_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/visits/data.stan"
  #include "modules/tumor/hyperparams.stan"
  #include "modules/state_space/data.stan"
  #include "modules/multistate/flags.stan"
  #include "modules/multistate/data.stan"
  #include "modules/multistate/hyperparams.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  int<lower = 0, upper = 1> fit_multistate_data;

  #include "modules/state_space/lfo_data.stan"
}

transformed data {
  #include "_hierarchy_transformed_data.stan"
  #include "_hmc_routing_transformed_data.stan"
  int max_all_t = max(max(t_patient_visits) + 1, extend_max_all_t);
  int<lower=0> max_t_width = max_all_t - min(t_patient_visits) + 1;
  #include "_visit_transformed_data.stan"
  #include "_full_model_transformed_data.stan"
  #include "modules/visits/transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"
  #include "modules/multistate/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/tumor/parameters.stan"
  #include "modules/multistate/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "modules/state_space/transformed_parameters.stan"
  #include "_ms_time_varying_covar.stan"
  #include "modules/multistate/transformed_parameters.stan"
}

generated quantities {
  #include "_lfo_endpoints_generated_quantities.stan"  
}