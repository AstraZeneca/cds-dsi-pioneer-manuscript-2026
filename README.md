```bash
quarto/
├── deck.qmd                              # revealjs slides deck. Not longer using this, but it still has some useful plots.
└── early-predict-bc-survival.qmd         # This is the primary notebook in which all analysis is done.
stan/
├── base_data.stan                        # This is the data section to be shared by all models.
├── pfs_functions.stan                    # Stan functions to be used in the PFS model.
├── pfs.stan                              # The main PFS model file.
├── tumor_data.stan                       # This is the tumor model settings data, useful for jointly fitting PFS and tumor models.
├── tumor_model.stan                      # The code in the model section of the tumor model.
├── tumor_parameters.stan                 # Parameters used in the tumor model.
├── tumor.stan                            # Main tumor model file (standalone).
├── tumor_transformed_data.stan           # Tranformed data section part for the tumor model.
├── tumor_transformed_parameters.stan     # Transformed parameters section for the tumor model.
└── util.stan                             # Utility functions to be shared by all models.
r/
├── download-entimice-data.R              # Script to download data.
├── sbc.R                                 # Standalone script to run simulation-based calibration.
├── sbc.sh                                # SLURM script to run sbc.R.
└── util.R                                # Shared functions.
```
