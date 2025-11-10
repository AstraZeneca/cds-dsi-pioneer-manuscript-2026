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
  #include "modules/other_events/data.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/other_events/hyperparams.stan"
  #include "modules/other_events/flags.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  #include "_sf-ssls-lfo-data.stan"
} 

transformed data {
  #include "../base_transformed_data.stan"
  #include "../tumor/tumor_transformed_data.stan"
  #include "_sf_transformed_data.stan"
  #include "modules/other_events/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/other_events/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
  #include "legacy/sf-ssls-parameters.stan"
}

transformed parameters {
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "_sf_transformed_parameters.stan"
  #include "modules/other_events/transformed_parameters.stan"
}

generated quantities {
  #include "_lfo_endpoints_generated_quantities.stan"  
}