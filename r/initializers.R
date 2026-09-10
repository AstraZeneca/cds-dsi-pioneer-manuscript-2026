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
    # Draw from inv_gamma prior (keeps mass away from zero)
    measure_sd_sld <- invgamma::rinvgamma(1, stan_data$measure_sd_sld_alpha, stan_data$measure_sd_sld_beta)

    # Hierarchical SDs (new naming convention)
    # Truncate at 0.05 to avoid near-zero inits that cause numerical issues
    tr_sd_patient_intercept <- pmax(0.05, abs(rnorm(1, 0, stan_data$tr_sd_patient_intercept_sd)))
    frac_sd_patient_intercept <- pmax(0.05, abs(rnorm(1, 0, stan_data$frac_sd_patient_intercept_sd)))

    # Existing initial state proportion parameters (new names)
    init_logit_loc_pop <- rnorm(1,
      stan_data$init_logit_loc_pop_mean,
      stan_data$init_logit_loc_pop_sd
    )
    init_sd_patient_intercept <- pmax(0.05, abs(rnorm(1, 0, stan_data$init_sd_patient_intercept_sd)))

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
      measure_sd_sld = measure_sd_sld
    )

    if (use_cross_process_corr) init_vals$L_process_corr <- L_process_corr

    init_vals$init_logit_loc_pop <- init_logit_loc_pop
    init_vals$init_sd_patient_intercept <- init_sd_patient_intercept
    init_vals$tr_sd_patient_intercept <- tr_sd_patient_intercept
    init_vals$frac_sd_patient_intercept <- frac_sd_patient_intercept

    if (isTRUE(stan_data$enable_static_init == 1L)) {
      init_vals$init_logit_static_loc_pop <- as.array(rnorm(
        1,
        stan_data$init_logit_static_loc_pop_mean,
        stan_data$init_logit_static_loc_pop_sd
      ))
    }

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
    pattern <- str_c("^", param_base, "\\[")
    matching_cols <- param_names %>%
      stringr::str_subset(pattern)

    # If no matches found, return NULL
    if (length(matching_cols) == 0) return(NULL)

    # Try to build the matrix
    result <- matrix(0, nrow = rows, ncol = cols)

    for (i in 1:rows) {
      for (j in 1:cols) {
        param <- str_c(param_base, "[", i, ",", j, "]")
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
    pattern <- str_c("^", param_base, "\\[")
    matching_cols <- param_names %>%
      stringr::str_subset(pattern)

    # If no matches found, return NULL
    if (length(matching_cols) == 0) return(NULL)

    # Try to build the vector
    result <- rep(NA, length)

    for (i in 1:length) {
      param <- str_c(param_base, "[", i, "]")
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
    "pop_log_growth_transition_rate", "measure_sd_sld",
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


#' Create SCLC Stan model initializer
#'
#' Returns function(chain_id) that draws starting values from priors.
#' Composes biomarker-specific (SLD tumor) initialization with the shared
#' multistate GP initialization from ms_init_values().
#'
#' @param stan_data Complete Stan data list (base data + priors + flags)
#' @return Function(chain_id) returning a named list of initial values
create_tumor_ssls_initializer <- function(stan_data) {
  .ms_init <- ms_init_values
  function(chain_id) {
    min_all_t <- min(stan_data$t_patient_visits)
    max_all_t <- max(max(stan_data$t_patient_visits) + 1, stan_data$extend_max_all_t %||% 0)
    max_t_width <- max_all_t - min_all_t + 1

    with(stan_data, {
      n_fgpl <- n_groups_per_level
      if (!is.null(n_forecast_patients)) n_fgpl[n_levels] <- n_forecast_patients

      n_raw_groups_tr_intercept  <- sum(n_fgpl[enable_level_intercept_tr %in% c(1L, 2L, 3L)])
      n_cp_groups_tr_intercept   <- sum(n_fgpl[enable_level_intercept_tr == 4L])
      n_raw_groups_frac_intercept <- sum(n_fgpl[enable_level_intercept_frac %in% c(1L, 2L, 3L)])
      n_cp_groups_frac_intercept  <- sum(n_fgpl[enable_level_intercept_frac == 4L])
      n_raw_groups_init_intercept <- sum(n_fgpl[enable_level_intercept_init %in% c(1L, 2L, 3L)])
      n_cp_groups_init_intercept  <- sum(n_fgpl[enable_level_intercept_init == 4L])

      n_raw_groups_tr_slope  <- sum(n_fgpl[(enable_level_intercept_tr %in% c(1L, 2L, 3L)) & (enable_level_cov_tr == 1L)])
      n_cp_groups_tr_slope   <- sum(n_fgpl[(enable_level_intercept_tr == 4L) & (enable_level_cov_tr == 1L)])
      n_raw_groups_frac_slope <- sum(n_fgpl[(enable_level_intercept_frac %in% c(1L, 2L, 3L)) & (enable_level_cov_frac == 1L)])
      n_cp_groups_frac_slope  <- sum(n_fgpl[(enable_level_intercept_frac == 4L) & (enable_level_cov_frac == 1L)])
      n_raw_groups_init_slope <- sum(n_fgpl[(enable_level_intercept_init %in% c(1L, 2L, 3L)) & (enable_level_cov_init == 1L)])
      n_cp_groups_init_slope  <- sum(n_fgpl[(enable_level_intercept_init == 4L) & (enable_level_cov_init == 1L)])

      rtruncnorm_actual <- function(n, sd, max_dev = 1.5) {
        raw <- rnorm(n, sd = sd * 0.6)
        pmax(-max_dev * sd, pmin(max_dev * sd, raw))
      }

      tr_sd_level <- pmax(0.05, abs(rnorm(n_levels, sd = tr_sd_level_intercept_sd)))
      frac_sd_level <- pmax(0.05, abs(rnorm(n_levels, sd = frac_sd_level_intercept_sd)))
      init_sd_level <- pmax(0.05, abs(rnorm(n_levels, sd = init_sd_level_intercept_sd)))

      tr_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_tr[lv] %in% c(1L, 2L, 3L)) rtruncnorm_actual(n_fgpl[lv], tr_sd_level[lv]) else NULL
      }))
      frac_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_frac[lv] %in% c(1L, 2L, 3L)) rtruncnorm_actual(n_fgpl[lv], frac_sd_level[lv]) else NULL
      }))
      init_level_dev <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (enable_level_intercept_init[lv] %in% c(1L, 2L, 3L)) rtruncnorm_actual(n_fgpl[lv], init_sd_level[lv]) else NULL
      }))

      raw_mask_tr   <- as.integer(enable_level_intercept_tr %in% c(1L, 2L, 3L))
      raw_mask_frac  <- as.integer(enable_level_intercept_frac %in% c(1L, 2L, 3L))
      raw_mask_init  <- as.integer(enable_level_intercept_init %in% c(1L, 2L, 3L))

      enabled_level_pos_tr <- c(1L, cumsum(n_fgpl * raw_mask_tr) + 1L)
      tr_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!raw_mask_tr[lv]) return(NULL)
        idx <- enabled_level_pos_tr[lv]:(enabled_level_pos_tr[lv + 1] - 1)
        tr_level_dev[idx] / tr_sd_level[lv]
      }))

      enabled_level_pos_frac <- c(1L, cumsum(n_fgpl * raw_mask_frac) + 1L)
      frac_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!raw_mask_frac[lv]) return(NULL)
        idx <- enabled_level_pos_frac[lv]:(enabled_level_pos_frac[lv + 1] - 1)
        frac_level_dev[idx] / frac_sd_level[lv]
      }))

      enabled_level_pos_init <- c(1L, cumsum(n_fgpl * raw_mask_init) + 1L)
      init_raw_level <- unlist(lapply(seq_len(n_levels), function(lv) {
        if (!raw_mask_init[lv]) return(NULL)
        idx <- enabled_level_pos_init[lv]:(enabled_level_pos_init[lv + 1] - 1)
        init_level_dev[idx] / init_sd_level[lv]
      }))

      tr_cp_level_intercept   <- rnorm(n_cp_groups_tr_intercept,   0, tr_sd_level[enable_level_intercept_tr   == 4L] * 0.2)
      frac_cp_level_intercept <- rnorm(n_cp_groups_frac_intercept, 0, frac_sd_level[enable_level_intercept_frac == 4L] * 0.2)
      init_cp_level_intercept <- rnorm(n_cp_groups_init_intercept, 0, init_sd_level[enable_level_intercept_init == 4L] * 0.2)

      tr_sd_level_slope <- lapply(seq_len(n_levels), function(lv) {
        pmax(0.01, abs(rnorm(n_covar, sd = tr_sd_level_slope_sd[[lv]])))
      })
      frac_sd_level_slope <- lapply(seq_len(n_levels), function(lv) {
        pmax(0.01, abs(rnorm(n_covar, sd = frac_sd_level_slope_sd[[lv]])))
      })
      init_sd_level_slope <- lapply(seq_len(n_levels), function(lv) {
        pmax(0.01, abs(rnorm(n_covar, sd = init_sd_level_slope_sd[[lv]])))
      })

      tr_raw_level_slope   <- matrix(rnorm(n_raw_groups_tr_slope   * n_covar, sd = 0.5),
                                     nrow = n_raw_groups_tr_slope,   ncol = n_covar)
      frac_raw_level_slope <- matrix(rnorm(n_raw_groups_frac_slope * n_covar, sd = 0.5),
                                     nrow = n_raw_groups_frac_slope, ncol = n_covar)
      init_raw_level_slope <- matrix(rnorm(n_raw_groups_init_slope * n_covar, sd = 0.5),
                                     nrow = n_raw_groups_init_slope, ncol = n_covar)

      tr_cp_level_slope   <- matrix(rnorm(n_cp_groups_tr_slope   * n_covar, sd = 0.1),
                                    nrow = n_cp_groups_tr_slope,   ncol = n_covar)
      frac_cp_level_slope <- matrix(rnorm(n_cp_groups_frac_slope * n_covar, sd = 0.1),
                                    nrow = n_cp_groups_frac_slope, ncol = n_covar)
      init_cp_level_slope <- matrix(rnorm(n_cp_groups_init_slope * n_covar, sd = 0.1),
                                    nrow = n_cp_groups_init_slope, ncol = n_covar)

      biomarker_init <- lst(
        tr_loc_pop = rnorm(1, tr_loc_pop_mean, tr_loc_pop_sd),
        frac_logit_loc_pop = rnorm(1, frac_logit_loc_pop_mean, frac_logit_loc_pop_sd),
        init_logit_loc_pop = rnorm(1, init_logit_loc_pop_mean, init_logit_loc_pop_sd),

        tr_sd_level_intercept_raw = as.array(tr_sd_level[enable_level_intercept_tr %in% c(2L, 4L)]),
        tr_raw_level_intercept = as.array(tr_raw_level),
        tr_cp_level_intercept = as.array(tr_cp_level_intercept),
        frac_sd_level_intercept_raw = as.array(frac_sd_level[enable_level_intercept_frac %in% c(2L, 4L)]),
        frac_raw_level_intercept = as.array(frac_raw_level),
        frac_cp_level_intercept = as.array(frac_cp_level_intercept),
        init_sd_level_intercept_raw = as.array(init_sd_level[enable_level_intercept_init %in% c(2L, 4L)]),
        init_raw_level_intercept = as.array(init_raw_level),
        init_cp_level_intercept = as.array(init_cp_level_intercept),

        tr_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_tr) array(rnorm(n_covar, tr_coef_qr_pop_mean, tr_coef_qr_pop_sd), dim = n_covar),
        frac_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_frac) array(rnorm(n_covar, frac_coef_qr_pop_mean, frac_coef_qr_pop_sd), dim = n_covar),
        init_coef_qr_pop = if (n_covar > 0 && enable_pop_cov_init) array(rnorm(n_covar, init_coef_qr_pop_mean, init_coef_qr_pop_sd), dim = n_covar),

        tr_sd_level_slope = if (n_covar > 0) tr_sd_level_slope,
        tr_raw_level_slope = if (n_covar > 0) tr_raw_level_slope,
        tr_cp_level_slope = if (n_covar > 0) tr_cp_level_slope,
        frac_sd_level_slope = if (n_covar > 0) frac_sd_level_slope,
        frac_raw_level_slope = if (n_covar > 0) frac_raw_level_slope,
        frac_cp_level_slope = if (n_covar > 0) frac_cp_level_slope,
        init_sd_level_slope = if (n_covar > 0) init_sd_level_slope,
        init_raw_level_slope = if (n_covar > 0) init_raw_level_slope,
        init_cp_level_slope = if (n_covar > 0) init_cp_level_slope,

        tr_raw_patient_process_noise = if (enable_patient_process_noise_tr) {
          matrix(0, nrow = n_patients, ncol = max_t_width)
        },
        tr_log_sd_pop_process_noise = if (enable_patient_process_noise_tr) array(rnorm(1, mean = log(0.05), sd = 0.5)),
        tr_sd_patient_log_sd_process_noise = if (enable_patient_process_noise_tr) array(abs(rnorm(1, sd = 0.3))),
        tr_raw_patient_log_sd_process_noise = if (enable_patient_process_noise_sd_tr) rep(0, n_patients),
        tr_logit_phi_pop_process_noise = if (enable_patient_process_noise_tr) array(rnorm(1, mean = 2, sd = 1)),
        tr_sd_patient_phi_process_noise = if (enable_patient_process_noise_tr) array(abs(rnorm(1, sd = 0.1))),
        tr_raw_patient_phi_process_noise = if (enable_patient_process_noise_phi_tr) rep(0, n_patients),

        tr_raw_pop_process_noise = if (enable_pop_process_noise_tr) {
          rnorm(max_t_width, 0, 0.1)
        },
        tr_log_sd_pop_process_noise_pop = if (enable_pop_process_noise_tr) array(rnorm(1, mean = log(0.05), sd = 0.5)),
        tr_logit_phi_pop_process_noise_pop = if (enable_pop_process_noise_tr) array(rnorm(1, mean = 1.4, sd = 0.3)),

        measure_sd_sld = invgamma::rinvgamma(1, measure_sd_sld_alpha, measure_sd_sld_beta),
      )

      c(biomarker_init, .ms_init(environment()))
    }) |> compact()
  }
}


# nolint end: object_usage_linter
