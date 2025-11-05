functions {
  #include "../util.stan"
  #include "../pos.stan"
  #include "../gp.stan"
  #include "../pfs_functions.stan"
  #include "../lfo.stan"
  #include "_sf_functions.stan"
  #include "../recist.stanfunctions"
}  

data {
  #include "../base_data.stan"
  #include "../tumor/base_data.stan"
  #include "_sf_outcomes_info.stan"

  #include "legacy/sf-ssls-hyperparam.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  #include "_sf-ssls-lfo-data.stan"
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "_sf_transformed_data.stan"
//   #include "other_events_transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
//   #include "other_events_parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
  #include "legacy/sf-ssls-parameters.stan"
}

transformed parameters {
//   #include "other_events_transformed_parameters.stan"
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "legacy/sf-ssls-transformed_parameters.stan"
}

generated quantities {
  #include "_lfo_endpoints_generated_quantities.stan"  
}