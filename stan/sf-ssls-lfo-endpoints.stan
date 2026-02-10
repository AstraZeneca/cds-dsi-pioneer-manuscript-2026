functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "state_space/functions.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "base_data.stan"
  #include "modules/tumor/data.stan"
  #include "state_space/data.stan"

  #include "modules/measurement/hyperparams.stan"
  #include "modules/other_events/data.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/other_events/hyperparams.stan"
  #include "modules/measurement/flags.stan"
  #include "modules/other_events/flags.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  #include "state_space/lfo_data.stan"
} 

transformed data {
  #include "base_transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/measurement/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "state_space/transformed_data.stan"
  #include "modules/other_events/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/measurement/parameters.stan"
  #include "modules/other_events/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/measurement/transformed_parameters.stan"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "state_space/transformed_parameters.stan"
  #include "modules/other_events/transformed_parameters.stan"
}

generated quantities {
  #include "_lfo_endpoints_generated_quantities.stan"  
}