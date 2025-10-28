# nolint start: object_usage_linter

# Stan Initializer Function Factories for CRCR and PFS Models
#
# This file contains functions that create initializer functions for Stan models,
# specifically for Competing Risks Confirmed Response (CRCR) and Progression-Free
# Survival (PFS) models. These initializers are designed to provide starting values
# for Stan MCMC sampling drawn from the prior distributions of the parameters,
# which can improve convergence and sampling efficiency.
#
# The main functions in this file are:
#
# 1. create_crcr_initializer:
#    Creates an initializer function for CRCR models. It generates initial values
#    from the prior distributions for parameters related to the baseline hazard,
#    including GP components.
#
# 2. create_crcr_pfs_initializer:
#    Creates an initializer function for combined CRCR and PFS models. It extends
#    the CRCR initializer with additional PFS-specific parameters, drawing initial
#    values from their respective prior distributions. This includes priors for
#    tumor stimulation effects and covariate effects.
#
# These functions take Stan data (which includes prior specifications) as input
# and return a function that generates initial values for a given chain. They
# account for various model configurations, such as separate baseline hazards,
# proportional hazards, and trial-level effects.

create_crcr_initializer <- function(stan_data, n_causes = 2) {
  base_sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  prop_sep_trial <- if (stan_data$separate_prop_hazard) stan_data$n_trials else 1 
  max_confresp_week <- max(max(stan_data$confirmed_response_week), stan_data$extend_max_confresp_week)
  
  function(chain_id) {
    init_vals <- lst(
      # array[n_causes] vector[n_base_separate_trials] log_crcr_lambda_gp_intercept;
      log_crcr_lambda_gp_intercept = matrix(
        with(stan_data, rnorm(n_causes * base_sep_trial, t(log_crcr_lambda_gp_intercept_mean), t(log_crcr_lambda_gp_intercept_sd))),
        ncol = base_sep_trial, byrow = FALSE 
      ),
      # array[n_causes] vector<lower = 0>[n_base_separate_trials] log_crcr_lambda_gp_alpha;
      log_crcr_lambda_gp_alpha = matrix(abs(rnorm(n_causes * base_sep_trial, sd = stan_data$log_crcr_lambda_gp_alpha_sd)), nrow = n_causes),

      # array[n_causes] vector<lower = 0>[n_base_separate_trials] log_crcr_lambda_gp_rho;
      log_crcr_lambda_gp_rho = with(stan_data, matrix(invgamma::rinvgamma(n_causes * base_sep_trial, log_crcr_lambda_gp_rho_alpha, log_crcr_lambda_gp_rho_beta), nrow = n_causes)),

      # array[n_causes] matrix[n_base_separate_trials, max_confresp_week] log_crcr_lambda_gp_eta;
      log_crcr_lambda_gp_eta = rnorm(n_causes * base_sep_trial * max_confresp_week) |> array(c(n_causes, base_sep_trial, max_confresp_week))
    )
    
    # if (!stan_data$no_prop_hazard) {
    #   init_vals <- init_vals |> 
    #     list_assign(
    #       # array[n_causes, n_prop_separate_trials] vector[no_prop_hazard ? 0 : n_tumor_covar] crcr_tumor_stim_pop_coef;
    #       crcr_tumor_stim_pop_coef = with(
    #         stan_data, 
    #         rnorm(n_causes * n_tumor_covar, 0, crcr_tumor_stim_pop_coef_sd) |>
    #           array(c(n_causes, prop_sep_trial, n_tumor_covar))  
    #       ),
    #       
    #       # array[n_causes, n_prop_separate_trials] vector[no_prop_hazard ? 0 : n_covar] crcr_covar_effect;  
    #       crcr_covar_coef = with(
    #         stan_data, 
    #         rnorm(n_causes * prop_sep_trial * n_covar, crcr_covar_effect_mean, crcr_covar_effect_sd) |>
    #           array(c(n_causes, prop_sep_trial, n_covar))  
    #       ),
    #     )
    #   
    #   if (stan_data$add_trial_level_prop_hazard) {
    #     init_vals <- init_vals |> 
    #       list_assign(
    #         # vector<lower = 0>[add_trial_level_prop_hazard && !no_prop_hazard ? n_tumor_covar + n_covar : 0] crcr_covar_trial_sd;
    #         crcr_covar_trial_sd = with(stan_data, rnorm(n_tumor_covar + n_covar, 0, crcr_covar_trial_sd_sd)),
    #         
    #         # array[n_causes, add_trial_level_prop_hazard && !no_prop_hazard ? n_trials : 0] vector[n_tumor_covar + n_covar] raw_crcr_covar_trial_coef;
    #         raw_crcr_covar_trial_coef = with(stan_data, rnorm(n_causes * prop_sep_trial * (n_tumor_covar + n_covar)) |> array(c(n_causes, prop_sep_trial, n_tumor_covar + n_covar)))
    #       )
    #   }
    #     
    # }
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |> 
        list_assign(
          # vector[add_trial_level_baseline_hazard ? n_causes : 0] log_crcr_lambda_gp_trial_intercept_sd;
          log_crcr_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(n_causes, sd = log_crcr_lambda_gp_trial_intercept_sd_sd)))
        )
    }
    
    return(init_vals) 
  }
}

create_crcr_pfs_initializer <- function(stan_data, n_causes = 2) {
  prop_hazard <- !stan_data$no_prop_hazard 
  crcr_init_fun <- create_crcr_initializer(stan_data, n_causes)
  baseline_sep_trial <- if (stan_data$separate_baseline_hazard) stan_data$n_trials else 1 
  prop_sep_trial <- if (stan_data$separate_prop_hazard) stan_data$n_trials else 1 
  max_all_t <- max(max(stan_data$t_measure) + 1, stan_data$extend_max_all_t)
  
  function(chain_id) {
    init_vals <- crcr_init_fun(chain_id) |> 
      list_assign(
        # vector<lower = 0>[n_base_separate_trials] log_lambda_gp_alpha;
        log_lambda_gp_alpha = abs(rnorm(baseline_sep_trial, stan_data$log_lambda_gp_alpha_sd)),
        
        # vector<lower = 0>[n_base_separate_trials] log_lambda_gp_rho;
        log_lambda_gp_rho = with(stan_data, invgamma::rinvgamma(baseline_sep_trial, log_lambda_gp_rho_alpha, log_lambda_gp_rho_beta)),
        
        # matrix[n_base_separate_trials, max_all_t] log_lambda_gp_eta;
        log_lambda_gp_eta = matrix(rnorm(baseline_sep_trial * max_all_t), nrow = baseline_sep_trial),
        
        # vector[n_base_separate_trials] log_lambda_gp_intercept;
        log_lambda_gp_intercept = with(stan_data, rnorm(baseline_sep_trial, log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd)),
        
        tumor_stim_pop_coef = if (prop_hazard) { 
          with(stan_data, matrix(rnorm(prop_sep_trial * n_tumor_covar, sd = c(t(tumor_stim_pop_coef_sd))), byrow = TRUE, nrow = prop_sep_trial))
        } else {
          array(dim = c(0, stan_data$n_tumor_covar))
        },
        covar_effect = if (prop_hazard) { 
          with(stan_data, matrix(rnorm(prop_sep_trial * n_covar, c(t(covar_effect_mean)), c(t(covar_effect_sd))), byrow = TRUE, nrow = prop_sep_trial))
        },
        conf_resp_effect = if (prop_hazard) { 
          with(stan_data, rnorm(prop_sep_trial, conf_resp_effect_mean, conf_resp_effect_sd))
        }
      ) |> 
      compact()
    
    if (with(stan_data, add_trial_level_baseline_hazard && !separate_baseline_hazard)) {
      init_vals <- init_vals |>  
        list_assign(
          log_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_intercept_sd_sd))),
          log_lambda_gp_trial_alpha = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_alpha_sd)))
        )  
    }
    
    if (with(stan_data, add_trial_level_prop_hazard && !separate_prop_hazard) && prop_hazard) {
      init_vals <- init_vals |>  
        list_assign(covar_trial_sd = with(stan_data, abs(rnorm(n_tumor_covar + n_covar + 1, sd = covar_trial_sd_sd))))  
    }
    
    return(init_vals)
  }
}

create_tumor_initializer <- function(stan_data) {
  max_all_t <- max(max(stan_data$t_measure) + 1, stan_data$extend_max_all_t)
  function(chain_id) {
    init_vals <- lst(
      pop_tumor_gp_alpha = abs(rnorm(1, sd = stan_data$pop_tumor_gp_alpha_sd)),
      
      log_pop_tumor_gp_rho = with(stan_data, rlnorm(1, pop_tumor_gp_rho_meanlog, pop_tumor_gp_rho_sdlog)),
      
      patient_tumor_intercept_sd = abs(rnorm(1, sd = stan_data$patient_tumor_intercept_sd_sd)), 
      raw_patient_tumor_intercept_effect = rnorm(stan_data$n_patients), 
      patient_tumor_intercept_effect = with(stan_data, rnorm(n_patients, sd = patient_tumor_intercept_sd)), 
    ) |> 
      compact()
    
    return(init_vals)
  } 
}

create_tumor_ss_initializer <- function(stan_data) {
  function(chain_id) {
    # Calculate number of visits minus 1 for training patients only
    n_total_train_visits_m1 <- sum(stan_data$n_patient_visits) - n_patients
    
    # Population-level parameters (using new naming convention)
    tr_loc_pop <- rnorm(1, stan_data$tr_loc_pop_mean, stan_data$tr_loc_pop_sd)
    frac_logit_loc_pop <- rnorm(1, stan_data$frac_logit_loc_pop_mean, stan_data$frac_logit_loc_pop_sd)

    # GP parameters
    log_pop_tumor_gp_rho <- rnorm(1, stan_data$pop_tumor_gp_rho_meanlog, stan_data$pop_tumor_gp_rho_sdlog)
    log_patient_tumor_gp_rho_sd <- abs(rnorm(1, 0, stan_data$log_patient_tumor_gp_rho_sd_sd))

    # Growth lag parameters
    pop_log_growth_lag <- rnorm(1, stan_data$growth_lag_mean, stan_data$growth_lag_sd)
    pop_log_growth_transition_rate <- abs(rnorm(1, 0, stan_data$log_growth_transition_rate_sd))
    patient_log_growth_lag_sd <- abs(rnorm(1, 0, stan_data$patient_log_growth_lag_sd_sd))

    # Process noise
    pop_process_sd <- c(
      abs(rnorm(1, 0, stan_data$pop_decrease_process_sd_sd)),
      abs(rnorm(1, 0, stan_data$pop_growth_process_sd_sd))
    )
    measure_sd <- abs(rnorm(1, 0, stan_data$measure_sd_sd))

    # Hierarchical SDs (new naming convention)
    tr_sd_patient_intercept <- abs(rnorm(1, 0, stan_data$tr_sd_patient_intercept_sd))
    frac_sd_patient_intercept <- abs(rnorm(1, 0, stan_data$frac_sd_patient_intercept_sd))

    # Existing initial state proportion parameters (new names)
    init_logit_loc_pop <- rnorm(1,
      stan_data$init_logit_loc_pop_mean,
      stan_data$init_logit_loc_pop_sd
    )
    init_sd_patient_intercept <- abs(rnorm(1, 0, stan_data$init_sd_patient_intercept_sd))

    use_cross_process_corr <- !stan_data$independ_cross_process_noise
    L_process_corr <- if (use_cross_process_corr) diag(2) else matrix(numeric(0), 0, 0)

    init_vals <- list(
      tr_loc_pop = tr_loc_pop,
      frac_logit_loc_pop = frac_logit_loc_pop,
      log_pop_tumor_gp_rho = log_pop_tumor_gp_rho,
      log_patient_tumor_gp_rho_sd = log_patient_tumor_gp_rho_sd,
      pop_log_growth_lag = pop_log_growth_lag,
      pop_log_growth_transition_rate = pop_log_growth_transition_rate,
      patient_log_growth_lag_sd = patient_log_growth_lag_sd,
      pop_process_sd = pop_process_sd,
      measure_sd = measure_sd
    )

    if (use_cross_process_corr) init_vals$L_process_corr <- L_process_corr

    init_vals$init_logit_loc_pop <- init_logit_loc_pop
    init_vals$init_sd_patient_intercept <- init_sd_patient_intercept
    init_vals$tr_sd_patient_intercept <- tr_sd_patient_intercept
    init_vals$frac_sd_patient_intercept <- frac_sd_patient_intercept

  if (!stan_data$enable_pop_cov_tr) {
      init_vals$tr_raw_patient_intercept <- rep(0, n_patients)
      init_vals$frac_raw_patient_intercept <- rep(0, n_patients)
    }

  if (!stan_data$enable_pop_cov_init) {
      init_vals$init_raw_patient_intercept <- rep(0, n_patients)
    }

  if (!stan_data$enable_pop_cov_tr) {
      init_vals$raw_patient_log_growth_lag <- rep(0, n_patients)
    }

  if (!stan_data$independ_long_process_noise && !stan_data$enable_pop_cov_tr) {
      init_vals$raw_log_patient_tumor_gp_rho_effect <- rep(0, n_patients)
    }

    init_vals$raw_patient_process_noise <- matrix(0, nrow = n_total_train_visits_m1, ncol = 2)
    init_vals$raw_states <- matrix(0, nrow = n_total_train_visits_m1, ncol = 2)
    init_vals
  }
}

create_tumor_ss_pathfinder_initializer <- function(pathfinder_fit, stan_data) {
  # Extract draws from the pathfinder fit
  draws_df <- posterior::as_draws_df(pathfinder_fit$draws())
  
  # Get parameter names
  param_names <- colnames(draws_df) |>  
    stringr::str_subset("^\\.", negate = TRUE) |>  
    stringr::str_subset("lp__|divergent__", negate = TRUE)
  
  # Calculate number of visits minus 1 for training patients only
  n_total_train_visits_m1 <- sum(stan_data$n_patient_visits) - n_patients
  
  # Flag for model configuration
  use_cross_process_corr <- !stan_data$independ_cross_process_noise
  use_long_process_corr <- !stan_data$independ_long_process_noise
  
  # Function to get parameter matrix from draws
  extract_matrix_param <- function(param_base, rows, cols) {
    pattern <- paste0("^", param_base, "\\[")
    matching_cols <- param_names %>% 
      stringr::str_subset(pattern)
    
    # If no matches found, return NULL
    if (length(matching_cols) == 0) return(NULL)
    
    # Try to build the matrix
    result <- matrix(0, nrow = rows, ncol = cols)
    
    for (i in 1:rows) {
      for (j in 1:cols) {
        param <- paste0(param_base, "[", i, ",", j, "]")
        if (param %in% param_names) {
          result[i, j] <- NA  # Just placeholder to check which elements exist
        }
      }
    }
    
    # Return NULL if empty matrix
    if (all(is.na(result))) return(NULL)
    
    return(result)
  }
  
  # Function to get parameter vector from draws
  extract_vector_param <- function(param_base, length) {
    pattern <- paste0("^", param_base, "\\[")
    matching_cols <- param_names %>% 
      stringr::str_subset(pattern)
    
    # If no matches found, return NULL
    if (length(matching_cols) == 0) return(NULL)
    
    # Try to build the vector
    result <- rep(NA, length)
    
    for (i in 1:length) {
      param <- paste0(param_base, "[", i, "]")
      if (param %in% param_names) {
        result[i] <- NA  # Just placeholder to check which elements exist
      }
    }
    
    # Return NULL if empty vector
    if (all(is.na(result))) return(NULL)
    
    return(result)
  }
  
  # Collect parameter information
  scalar_params <- c(
    "tr_loc_pop", "frac_logit_loc_pop",
    "log_pop_tumor_gp_rho", "pop_log_growth_lag",
    "pop_log_growth_transition_rate", "measure_sd",
    "init_logit_loc_pop",
    "log_patient_tumor_gp_rho_sd", "tr_sd_patient_intercept",
    "patient_log_growth_lag_sd", "init_sd_patient_intercept",
    "frac_sd_patient_intercept"
  ) |> purrr::keep(\(p) p %in% param_names)
  
  # Vector parameters
  vector_params <- list(
    pop_process_sd = 2
  ) 
  
  # Patient-level parameters
  patient_params <- list(
    tr_raw_patient_intercept = if (!stan_data$enable_pop_cov_tr) n_patients,
    frac_raw_patient_intercept = if (!stan_data$enable_pop_cov_tr) n_patients,
    raw_patient_log_growth_lag = if (!stan_data$enable_pop_cov_tr) n_patients,
    init_raw_patient_intercept = if (!stan_data$enable_pop_cov_init) n_patients,
    raw_log_patient_tumor_gp_rho_effect = if (use_long_process_corr && !stan_data$enable_pop_cov_tr) n_patients
  ) |> compact() 
  
  # Matrix parameters
  matrix_params <- list(
    L_process_corr = if (use_cross_process_corr) c(2, 2),
    raw_patient_process_noise = c(n_total_train_visits_m1, 2)
  ) |> 
    compact()
  
  # Return the initializer function
  function(chain_id) {
    # Randomly select a draw
    # draw_idx <- sample(1:nrow(draws_df), 1)
    draw <- draws_df |> sample_n(1)
    
    # Add scalar parameters
    init_vals <- scalar_params |>  
      map(\(p) as.numeric(pull(draw, p))) |> 
      set_names(scalar_params) 
    
    # Make sure pop_log_rate_ratio meets constraint if it exists
    if (!is.null(init_vals$pop_log_rate_ratio)) {
      init_vals$pop_log_rate_ratio <- max(init_vals$pop_log_rate_ratio, 0.125)
    }
    
    init_vals <- imap(c(vector_params, patient_params), \(s, p) unlist(draw[1, str_glue("{p}[{seq(s)}]")], use.names = FALSE)) |> 
      c(init_vals)
    
    init_vals <- imap(matrix_params, function(s, p) {
      param <- crossing(!!!map(s, seq)) |> 
        set_names(c("i", "j")) %$% 
        str_glue("{p}[{i},{j}]")
      
      unlist(draw[1, param]) |> 
        matrix(s[1], s[2], byrow = TRUE)
    }) |> 
      c(init_vals)
    
    # Default to simple initializers if not found in pathfinder results
    
    # Make sure process_sd is initialized if not already
    if (is.null(init_vals$pop_process_sd)) {
      init_vals$pop_process_sd <- c(
        abs(rnorm(1, 0, stan_data$pop_decrease_process_sd_sd)),
        abs(rnorm(1, 0, stan_data$pop_growth_process_sd_sd))
      )
    }
    
    # Process correlation matrix if needed
    if (use_cross_process_corr && is.null(init_vals$L_process_corr)) {
      init_vals$L_process_corr <- diag(2)
    }
    
    # Make sure raw states and process noise are initialized
    if (is.null(init_vals$raw_states)) {
      init_vals$raw_states <- matrix(0, nrow = n_total_train_visits_m1, ncol = 2)
    }
    
    if (is.null(init_vals$raw_patient_process_noise)) {
      init_vals$raw_patient_process_noise <- matrix(0, nrow = n_total_train_visits_m1, ncol = 2)
    }
    
    return(init_vals)
  }
}

create_tumor_ssls_initializer <- function(stan_data) {
  function(chain_id) {
    with(stan_data, {
      lst(
        tr_sd_trial_intercept = abs(rnorm(1, sd = tr_sd_trial_intercept_sd)),
        tr_raw_trial_intercept = if (enable_trial_intercept_tr) rnorm(n_trials),
        frac_sd_trial_intercept = abs(rnorm(1, sd = init_sd_trial_intercept_sd)),
        frac_raw_trial_intercept = if (enable_trial_intercept_frac) rnorm(n_trials),
        init_sd_trial_intercept = abs(rnorm(1, sd = init_sd_trial_intercept_sd)),
        init_raw_trial_intercept = if (enable_trial_intercept_init) rnorm(n_trials),

        # Patient-level hierarchical standard deviations (new naming)
        tr_sd_patient_intercept = abs(rnorm(1, sd = tr_sd_patient_intercept_sd)),
        tr_raw_patient_intercept = if (enable_patient_intercept_tr) rnorm(n_patients),
        frac_sd_patient_intercept = abs(rnorm(1, sd = frac_sd_patient_intercept_sd)),
        frac_raw_patient_intercept = if (enable_patient_intercept_frac) rnorm(n_patients),
        init_sd_patient_intercept = abs(rnorm(1, sd = init_sd_patient_intercept_sd)),
        init_raw_patient_intercept = if (enable_patient_intercept_init) rnorm(n_patients),

        tr_coef_qr_pop = if (n_covar > 0 && enable_trial_cov_tr) rnorm(n_covar, tr_coef_qr_pop_mean, tr_coef_qr_pop_sd),
        frac_coef_qr_pop = if (n_covar > 0 && enable_trial_cov_frac) rnorm(n_covar, frac_coef_qr_pop_mean, frac_coef_qr_pop_sd),
        init_coef_qr_pop = if (n_covar > 0 && enable_trial_cov_init) rnorm(n_covar, init_coef_qr_pop_mean, init_coef_qr_pop_sd),

        tr_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_tr) abs(rnorm(n_covar, sd = tr_sd_trial_slope_sd)),
        tr_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_tr) matrix(rnorm(n_covar * n_trials), n_trials, n_covar),
        frac_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_frac) abs(rnorm(n_covar, sd = frac_sd_trial_slope_sd)),
        frac_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_frac) matrix(rnorm(n_covar * n_trials), n_trials, n_covar),
        init_sd_trial_slope = if (n_covar > 0 && enable_trial_cov_init) abs(rnorm(n_covar, sd = init_sd_trial_slope_sd)),
        init_raw_trial_slope = if (n_covar > 0 && enable_trial_cov_init) matrix(rnorm(n_covar * n_trials), n_trials, n_covar),

        # Patient-level slope SDs and raw effects for tr, frac, and init modules
        tr_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_tr) abs(rnorm(n_covar, sd = tr_sd_patient_slope_sd)),
        tr_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_tr) matrix(rnorm(n_patients * n_covar), nrow = n_patients, ncol = n_covar),
        frac_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_frac) abs(rnorm(n_covar, sd = frac_sd_patient_slope_sd)),
        frac_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_frac) matrix(rnorm(n_patients * n_covar), nrow = n_patients, ncol = n_covar),
        init_sd_patient_slope = if (n_covar > 0 && enable_patient_cov_init) abs(rnorm(n_covar, sd = init_sd_patient_slope_sd)),
        init_raw_patient_slope = if (n_covar > 0 && enable_patient_cov_init) matrix(rnorm(n_patients * n_covar), nrow = n_patients, ncol = n_covar),
      )
    }) |> compact()
  }
}

# nolint end: object_usage_linter