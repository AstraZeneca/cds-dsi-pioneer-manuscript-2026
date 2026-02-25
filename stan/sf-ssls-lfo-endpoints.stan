functions {
  #include "util.stanfunctions"
  #include "pos.stanfunctions"
  #include "gp.stanfunctions"
  #include "pfs.stanfunctions"
  #include "lfo.stanfunctions"
  #include "modules/state_space/sf.stanfunctions"
  #include "modules/tumor/tumor.stanfunctions"
}

data {
  #include "_base_data.stan"
  #include "modules/tumor/data.stan"
  #include "modules/tumor/hyperparams.stan"
  #include "modules/state_space/data.stan"

  #include "modules/other_events/data.stan"
  #include "modules/tr/hyperparams.stan"
  #include "modules/frac/hyperparams.stan"
  #include "modules/init/hyperparams.stan"
  #include "modules/other_events/hyperparams.stan"
  #include "modules/other_events/flags.stan"
  #include "modules/tr/flags.stan"
  #include "modules/frac/flags.stan"
  #include "modules/init/flags.stan"

  #include "modules/state_space/lfo_data.stan"
} 

transformed data {
  #include "_base_transformed_data.stan"
  #include "modules/tumor/transformed_data.stan"
  #include "modules/tr/transformed_data.stan"
  #include "modules/frac/transformed_data.stan"
  #include "modules/init/transformed_data.stan"
  #include "modules/state_space/transformed_data.stan"

  // Generic biomarker baseline for other_events module (SLD mode)
  vector[n_patients] log_baseline_biomarker = log_baseline_sld;

  #include "modules/other_events/transformed_data.stan"
  #include "_lfo_transformed_data.stan"
}

parameters {
  #include "modules/tumor/parameters.stan"
  #include "modules/other_events/parameters.stan"
  #include "modules/tr/parameters.stan"
  #include "modules/frac/parameters.stan"
  #include "modules/init/parameters.stan"
}

transformed parameters {
  #include "modules/tr/transformed_parameters.stan"
  #include "modules/frac/transformed_parameters.stan"
  #include "modules/init/transformed_parameters.stan"
  #include "modules/state_space/transformed_parameters.stan"
  #include "modules/other_events/transformed_parameters.stan"
}

generated quantities {
  #include "_lfo_endpoints_generated_quantities.stan"  
}