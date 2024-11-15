This repository is used by PIONEER projects for predicting PFS and ORR outcomes in small and/or ongoing trials using data from older trials. Currently, there are two projects:

1. Predicting outcomes from one of the Breast-01 to -04 trials using the other three.
2. Predicting the LUNG trial outcomes from the Endometrial trial.

We predict outcomes using a Bayesian multilevel piece-wise constant proportional hazard survival model for PFS and, similarly, a competing risks model for confirmed response status (which we use as predictor in the PFS model).

# Repo Structure

``` bash
.
├── endometrial_to_lung_targets.R                                        # Endometrial to LUNG project targets file
├── _targets.R                                                        # Obsolete
└── quarto                                                            # Notebooks
    ├── confirmed-response-pfs.qmd                                    # Breast notebook
    ├── deck.qmd                                                      # Old deck
    ├── endometrial-to-lung.qmd                                         # Endometrial to LUNG notebook
    └── early-predict-bc-survival.qmd                                 # Old Breast notebook
├── r
│   ├── crcr.R                                                        # Common CRCR utility/posterior functions
│   ├── breast                                                # Breast project specific functions
│   │   ├── bg.R                             
│   │   ├── crcr.R
│   │   ├── posterior.R
│   │   ├── prepare_analysis_data.R
│   │   ├── sbc.R
│   │   └── util.R
│   ├── endometrial-to-lung                                             # Endometrial to LUNG specific functions (see similarly names files below)
│   │   ├── prepare_analysis_data.R
│   │   └── priors.R
│   ├── entimice                                                      # Scripts to download analysis datasets
│   │   ├── download-lung-entimice-data.R
│   │   ├── download-breast-entimice-data.R
│   │   ├── download-endometrial-entimice-data.R
│   │   └── entimice_functions.R
│   ├── initializers.R                                                # Common Stan initializers
│   ├── plot_functions.R                                              # Common plot functions
│   ├── posterior.R                                                   # Common functions for extracting samples from Stan fit objects
│   ├── prepare_analysis_data.R                                       # Common functions used for analysis and Stan data preparation
│   ├── priors.R                                                      # Common prior specification
│   ├── table_functions.R                                             # Common functions to generate {gt} tables
│   └── util.R                                                        # Common utility functions
└── stan                                                              # Stan statistical models folder
    ├── base_data.stan                                                # Data shared between all models
    ├── baseline_hazard                                               # GP PFS baseline hazard files
    │   ├── baseline_hazard_hyperparam.stan                           # Hyperparameters for GP priors
    │   ├── baseline_hazard_log_lik_prior_sense.stan                  # Calculations for {priorsense} and {loo}
    │   ├── baseline_hazard_parameters.stan                           
    │   ├── baseline_hazard_priors.stan                               
    │   └── baseline_hazard_transformed_parameters.stan               
    ├── base_transformed_data.stan                                    # Transformed data shared between all model
    ├── bootstrap                                                     # Bootstrap files
    │   ├── insample_bootstrap_data.stan                              
    │   ├── insample_bootstrap_gen_quants.stan                        
    │   ├── leave_out_trial_bootstrap_data.stan                       
    │   ├── leave_out_trial_bootstrap_functions.stan                  
    │   ├── leave_out_trial_bootstrap_gen_quants.stan                 
    │   └── leave_out_trial_bootstrap_transformed_data.stan           
    ├── crcr                                                          # Confirmed response competing risks (CRCR) model files
    │   ├── confresp-comprisk.stan                                    # Standalone CRCR model
    │   ├── crcr_baseline_hazard_hyperparam.stan                      
    │   ├── crcr_data.stan                                            
    │   ├── crcr_functions.stan
    │   ├── crcr_gen_quants.stan
    │   ├── crcr_hyperparam.stan
    │   ├── crcr_log_lik_prior_sense.stan                             # CRCR calculations needed for {priorsense} and {loo}
    │   ├── crcr_parameters.stan
    │   ├── crcr_priors.stan
    │   ├── crcr_transformed_data.stan
    │   └── crcr_transformed_parameters.stan
    ├── breast                                                # Breast project specific files
    │   ├── pfs2.stan
    │   ├── pfs_functions.stan
    │   ├── pfs_generated_quant.stan
    │   ├── pfs_orr.stan
    │   ├── pfs.stan
    │   ├── tumor_stim_hyperparam.stan
    │   ├── tumor_stim_parameters.stan
    │   ├── tumor_stim_priors.stan
    │   ├── tumor_stim_transformed_data.stan
    │   └── tumor_stim_transformed_parameters.stan
    ├── extern_pfs_functions.stan                                    # External PFS model functions (should be merged with pfs_functions.stan)
    ├── extern_util.stan                                             # External general utility function (should be merged with util.stan)
    ├── fixed_bootstrap_transformed_data.stan                       
    ├── pfs-confirmed-response.stan                                  # Main CRR + PFS model
    ├── pfs_functions.stan                                           # PFS specific functions
    ├── pfs_transformed_data.stan                            
    ├── recruit                                                      # Files for modelling trial recruitment rates and timing
    │   ├── recruit_parameters.stan
    │   ├── recruit_priors.stan
    │   ├── recruit_sample_maturity.stan
    │   └── recruit.stan                                             # Standalone model
    ├── tumor                                                        # Tumor specific modelling files
    │   ├── tumor_data.stan
    │   ├── tumor_model.stan                                         # Priors and likelihood calculations
    │   ├── tumor_parameters.stan
    │   ├── tumor.stan                                               # Standalone model
    │   └── tumor_transformed_parameters.stan
    └── util.stan                                                    # General utility functions
```

# Setup

1. Use `r/download-entimice-data.R` to download the studies data.
2. Render or execute the code in `quarto/early-predict-bc-survival.qmd`. This file is divided into two main sections:
    i. A tumor-size dynamics model using Gaussian processes.
    ii. The main survival model conditional on early tumor-size assessments.
3. There are a few supplementary plots also generated in `quarto/deck.qmd`.

# Projects

## Breast-01 to -04

## Predicting LUNG from Endometrial
