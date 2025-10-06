  #
  # This file contains functions for Leave-Future-Out (LFO) cross-validation,
  # model stacking, and various model evaluation metrics for survival analysis
  # and clinical trial data.}

# nolint start: object_usage_linter

get_trial_loo <- function(res, log_lik_var = "trial_log_lik", moment_match = TRUE, ...) {
  res$loo(log_lik_var, moment_match = moment_match, save_psis = TRUE, ...)
}

#' Add Stacked Results to Existing Results Data
#'
#' This function combines multiple model results using stacking weights to create
#' a new "stacked" model result. It adds this stacked result to the original data.
#'
#' @param res_data A data frame containing the original results data
#' @param stacking_weights A data frame containing the stacking weights for each model
#' @param ... Names of columns in res_data that contain the results to be stacked
#' @param by A character vector specifying the columns to join by (default: c("model_type", "trial"))
#'
#' @return A data frame containing the original results plus the stacked results
#'
#' @details
#' The function performs the following steps:
#' 1. Joins the original results with the stacking weights.
#' 2. Converts the weights to numeric values.
#' 3. For each combination of fit_type and trial, it creates a new "stacked" model result
#'    by combining the specified result columns using the stacking weights.
#' 4. Adds these stacked results to the original data.
#'
#' The stacking is performed using a nested function `stack_draws`, which:
#' - Allocates draws to each model based on its weight using `simplex_allocate`.
#' - Resamples the draws for each model according to its allocation.
#' - Combines the resampled draws into a single random variable (rvar).
#'
#' This approach ensures that the stacked result maintains the correct proportions
#' of draws from each model as specified by the stacking weights.
#'
add_stacked_results <- function(res_data, stacking_weights, ..., by = c("model_type", "trial")) {
  stack_draws <- function(rvs, simplex) {
    map2(rvs, simplex_allocate(simplex, ndraws(rvs[1])), \(rv, n) resample_draws(rv, ndraws = n)) |> 
      map(\(d) as.vector(draws_of(d))) |> 
      purrr::flatten_dbl() |> 
      rvar()
  }
  
  inner_join(res_data, stacking_weights, by = by) |>  
    mutate(weight = as.numeric(weight)) %>%
    bind_rows(
      filter(., !is.na(weight)) |> 
        group_by(fit_type, trial) |> 
        summarize(model_type = "stacked", across(c(...), \(res) stack_draws(res, weight))) 
    )
}

get_trial_c_index <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
  spread_rvars(res, trial_c_index[trial])
}

get_patient_pointwise_loo <- function(model_loo, stan_data) {
  as_tibble(stan_data[c("patient", "patient_trial")]) |> 
    rename(trial = patient_trial) |> 
    mutate(
      pareto_k_influence = loo::pareto_k_influence_values(model_loo), imputed = row_number() %in% stan_data$imputed_patients,
      elpd_loo = loo::pointwise(model_loo, "elpd_loo")
    )
}

#' Allocate Integer Values Based on a Simplex
#'
#' This function takes a simplex (a vector of non-negative values that sum to 1)
#' and a total integer value, and allocates the total across the simplex elements
#' as integer values, maintaining the proportions as closely as possible.
#'
#' @param simplex A numeric vector representing a simplex (sum should be 1)
#' @param total An integer representing the total to be allocated
#'
#' @return An integer vector of the same length as `simplex`, representing the allocation
#'
#' @details
#' The function works as follows:
#' 1. It first validates the inputs to ensure the simplex sums to 1 (within a small tolerance)
#'    and that the total is an integer.
#' 2. It then performs an initial allocation by multiplying each simplex value by the total
#'    and taking the floor of the result.
#' 3. If there's any remainder after the initial allocation, it distributes the remaining
#'    units one by one to the elements with the largest fractional parts.
#'
#' This approach ensures that the allocation is as close as possible to the proportions
#' specified by the simplex, while still resulting in integer values that sum to the specified total.
#'
#' @examples
#' simplex_allocate(c(0.3, 0.5, 0.2), 10)  # Returns c(3, 5, 2)
#' simplex_allocate(c(0.33, 0.33, 0.34), 100)  # Returns c(33, 33, 34)
#'
simplex_allocate <- function(simplex, total) {
  # Input validation
  if (abs(sum(simplex) - 1) > 1e-5) {
    stop("Input vector must sum to 1")
  }
  if (total %% 1 != 0) {
    stop("total must be an integer")
  }
  
  # Initial allocation using floor after multiplication
  raw_allocation <- simplex * total
  initial_allocation <- floor(raw_allocation)
  
  # Calculate remaining amount to distribute
  remainder <- total - sum(initial_allocation)
  
  if (remainder > 0) {
    # Get fractional parts
    fractional_parts <- raw_allocation - initial_allocation
    # Get indices that would sort in descending order
    sorted_indices <- order(fractional_parts, decreasing = TRUE)
    
    # Only distribute up to the remainder amount
    result <- initial_allocation
    if (remainder > 0) {
      result[sorted_indices[1:remainder]] <- result[sorted_indices[1:remainder]] + 1
    }
    return(result)
  } else {
    return(initial_allocation)
  }
}

get_loo_admin_brier_score <- function(res, loo_obj) {
  res |> 
    spread_draws(trial_admin_brier_score[i, t]) |> 
    ungroup() |> 
    select(.draw, i, t, trial_admin_brier_score) |> 
    pivot_wider(id_cols = c(.draw, t), names_from = i, values_from = trial_admin_brier_score) |> 
    select(!.draw) |> 
    nest(draws_matrix = !t) |> 
    transmute(t, mean_brier_score = map_dbl(draws_matrix, \(m) sum(loo::E_loo(as.matrix(m), loo_obj$psis_object, type = "mean")$value)))
}

# LFO functions ######

#' Generate Leave-Future-Out (LFO) Cutoff Days
#'
#' This function generates a sequence of cutoff dates and corresponding calendar days
#' for Leave-Future-Out (LFO) cross-validation in time series or longitudinal data analysis.
#'
#' @param first_cutoff_date Date. The starting date for the cutoff sequence.
#' @param last_date Date. The ending date for the cutoff sequence.
#' @param first_cutoff_day_idx Numeric. The calendar day index corresponding to the first cutoff date.
#' @param days_increment Numeric. The number of days between each cutoff.
#'
#' @return A tibble with columns:
#'   \item{n}{Integer. Sequential number for each cutoff.}
#'   \item{cutoff_date}{Date. The date of each cutoff.}
#'   \item{cutoff_calendar_day}{Numeric. The calendar day index for each cutoff.}
#'
#' @details
#' The function performs the following steps:
#' 1. Calculates the number of cutoffs based on the time span and increment.
#' 2. Generates a sequence of cutoff dates using the specified increment.
#' 3. Calculates the corresponding calendar day index for each cutoff date.
#'
#' This is particularly useful for setting up Leave-Future-Out cross-validation
#' in time-dependent analyses, such as clinical trials or longitudinal studies.
#'
get_lfo_cutoff_days <- function(first_cutoff_date, last_date, first_cutoff_day_idx, days_increment) {
  len <- time_length(last_date - first_cutoff_date, unit = "days") %/% days_increment + 1
  n <- seq(len) 
  
  tibble(
    n,
    cutoff_date = accumulate(n[-len], \(prev, n) prev + days(days_increment), .init = first_cutoff_date),
    cutoff_calendar_day = first_cutoff_day_idx + time_length(cutoff_date - first_cutoff_date, unit = "days"),
  ) 
}

#' Perform Leave-Future-Out (LFO) Cross-Validation
#'
#' This function implements Leave-Future-Out cross-validation for time series or longitudinal data,
#' particularly useful for Bayesian models fit with Stan.
#'
#' @param stan_data List. The data to be passed to the Stan model.
#' @param model Stan model object. The model to be fit.
#' @param cutoffs Data frame. Contains the cutoff points for LFO, typically generated by `get_lfo_cutoff_days`.
#' @param output_path String. Path where output files will be saved.
#' @param basename String. Base name for output files.
#' @param output_timestamp Logical. Whether to include a timestamp in output filenames (default: FALSE).
#' @param refit_n Integer. The initial cutoff point to start fitting from (default: min(cutoffs$n)).
#' @param k_threshold Numeric. Threshold for Pareto k statistic to determine when to refit (default: 0.7).
#' @param lean Logical. If TRUE, return only essential columns in the results (default: FALSE).
#' @param verbose Logical. If TRUE, print progress information (default: FALSE).
#' @param exact Logical. If TRUE, refit at every time point regardless of k statistic (default: FALSE).
#' @param fit_only Logical. If TRUE, only fit the model without computing LFO (default: FALSE).
#' @param iter_warmup Integer. Number of warmup iterations for Stan (default: 300).
#' @param iter_sampling Integer. Number of sampling iterations for Stan (default: 500).
#' @param ... Additional arguments passed to `sample_and_save`.
#'
#' @return A data frame containing LFO results, including log-likelihood and Pareto k statistics.
#'
#' @details
#' The function performs the following steps:
#' 1. Fits the model using the provided data and cutoffs.
#' 2. Computes LFO log-likelihood and Pareto k statistics.
#' 3. Determines if refitting is necessary based on the k statistic and threshold.
#' 4. Recursively calls itself with updated cutoffs if refitting is needed.
#'
#' This implementation allows for adaptive refitting, where the model is only refit
#' when the approximation quality (as measured by the Pareto k statistic) degrades.
#'
lfo <- function(
    stan_data, model, cutoffs, all_cutoffs, output_path, basename, initializer, 
    output_timestamp = FALSE, refit_n = min(cutoffs$n), 
    k_threshold = 0.7, lean = FALSE, verbose = FALSE, exact = FALSE, fit_only = FALSE, 
    iter_warmup = 300, iter_sampling = 500, parallel_chains = 4, adapt_delta = 0.9, future_window = 1, ...) {
  if (verbose) {
    cat("Starting on:\n")
    print(cutoffs)
    cat("\n")
  }
  
  remaining_cutoffs <- cutoffs |> filter(n >= refit_n) 
  remaining_all_cutoffs <- all_cutoffs |> filter(n >= refit_n)
  
  fit <- stan_data |>
    list_assign(cutoff_calendar_day = remaining_all_cutoffs$cutoff_calendar_day, n_cutoffs = nrow(remaining_all_cutoffs)) %>%
    sample_and_save(
      model,
      .,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling, parallel_chains = parallel_chains, adapt_delta = adapt_delta,
      init = initializer,
      output_dir = file.path(output_path, "fit"), output_basename = str_glue("{basename}-{refit_n}"),
      timestamp = output_timestamp, 
      ...
    ) 
  
  psis_results <- fit |> 
    lfo_log_lik(future_window = future_window) |> 
    mutate(across(c(n, m), \(x) x + refit_n - 1)) |> 
    left_join(select(remaining_all_cutoffs, n, cutoff_date, cutoff_calendar_day), by = "n") |> 
    mutate(tar_group = first(cutoffs$tar_group %||% NA_integer_), refit_n)
  
  if (fit_only) {
    return(lst(fit, psis_results))
  }
  
  if (lean) {
    psis_results <- psis_results |>
      select(n, m, contains("E_"))
  }
  
  next_cutoffs <- psis_results |> 
    filter(!is.na(k), k > k_threshold | exact, n > refit_n) %>%
    semi_join(remaining_cutoffs, ., by = "n")

  if (verbose) {
    cat("LFO results:\n")
    print(select(psis_results, n, m, refit_n, k))
    cat("\n")
  }
  
  if (nrow(next_cutoffs) > 0) {
    next_results <- lfo(
        stan_data, model, cutoffs, all_cutoffs, output_path, basename, initializer, output_timestamp, refit_n = min(next_cutoffs$n), 
        k_threshold, lean, verbose, exact, fit_only, iter_warmup, iter_sampling, parallel_chains, adapt_delta, ...
      )

    return(bind_rows(psis_results, next_results))
  } else {
    return(psis_results)
  }
}

lfo_drop_bad_approx <- function(lfo_res) {
  lfo_res |> 
    group_by(n) %>% 
    filter(if (has_name(., "refit_n")) min_rank(refit_n) == n() else TRUE) |> 
    ungroup()  
}

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

lfo_log_lik <- function(res, max_n = Inf, future_window = 1) {
  res |> 
    spread_rvars(patient_log_lik[n, m, i], patient_pfs_log_lik[n, m, i], patient_crcr_log_lik[n, m, i]) |>
    lfo_log_lik_rvar(max_n, future_window)
}

psis_resample <- function(l, w, recalc_full = FALSE) { #, negative_only = TRUE) {
  map2(l, w, function(ln, wn) {
    if (!is_null(ln)) {
      if (!is_null(wn)) {
        plyr::aaply(ln, 2, \(lni) log_sum_exp(lni + wn * all(wn < 0))) 
      } else if (recalc_full) { 
        plyr::aaply(ln, 2, log_mean_exp)
      }
    }
  })
}

#' Calculate Log-Likelihood for Leave-Future-Out (LFO) Cross-Validation
#'
#' This function processes log-likelihood values for Leave-Future-Out cross-validation,
#' computing various statistics including PSIS (Pareto Smoothed Importance Sampling) estimates.
#'
#' @param log_lik_rvar A data frame containing log-likelihood values as random variables (rvars),
#'        typically output from a Stan model. Expected columns include 'n', 'm', and various
#'        'patient_log_lik' columns.
#' @param max_n Integer. The maximum 'n' value to process (default: Inf).
#' @param future_window Integer. The number of future time points to consider (default: 1).
#'
#' @return A data frame with processed log-likelihood values and related statistics, including:
#'   \item{n}{Time index}
#'   \item{mean_patient_log_lik}{Mean log-likelihood}
#'   \item{psis_log_ratio}{PSIS estimates for log ratios}
#'   \item{k}{Pareto k values}
#'   \item{lwt}{PSIS weights}
#'   \item{approx_E_patient_log_lik}{Approximated expected log-likelihood}
#'   \item{E_patient_log_lik}{Expected log-likelihood}
#'   ... and similar columns for PFS and CRCR specific log-likelihoods
#'
#' @details
#' The function performs several steps:
#' 1. Filters and groups the input data.
#' 2. Calculates mean log-likelihoods and log ratios.
#' 3. Computes PSIS estimates and related statistics (k values, weights).
#' 4. Calculates approximated expected log-likelihoods using PSIS resampling.
#' 5. Renames and reorganizes columns for clarity.
#'
#' This function is crucial for assessing model performance in a time-series context.
#'
lfo_log_lik_rvar <- function(log_lik_rvar, max_n = Inf, future_window = 1) {
  log_lik_rvar |>   
    filter(m >= n) |>
    group_by(n, m) |> 
    summarize(across(matches("^patient(_.+)?_log_lik"), \(l) list(draws_of(l))), .groups = "drop") |>
    (function(d) {
      inner_join(
        filter(d, m == max(m)) |> select(!m), # From n to max(m), this is the out of sample loglik. For n = 1, that is the exact SAP.
        filter(d, n == 1) |> select(!n),      # From 1 to, this is the loglik for the additional periods of time that we want to PSIS to approximate.
                                              # This is relevant to predicting the _next_ row down.
        by = c("n" = "m"), suffix = c("", "_log_ratio")
      ) |>
        # This add loglik columns for M-SAP, rather than the full SAP we get from the above join.  
        left_join(
          # mutate(d, m = m - future_window + 1) |> filter(n == m), 
          filter(d, n == m - future_window + 1),
          by = "n", 
          suffix = c("", "_w") # _w is in reference to the m-sap "window"
        )
    })() |>
    filter(n <= max_n) |> 
    mutate(
      # fit = map(min_rank(n), \(nr) if (nr == 1) res),
      across(
        matches("^patient(_.+)?_log_lik(_w)?$"), 
        \(l) map_if(l, \(ln) !is_null(ln), \(ln) plyr::aaply(ln, 2, \(lni) log_mean_exp(lni))), 
        .names = "mean_{.col}"
      ),
      across(
        matches("^patient(_.+)?_log_lik_log_ratio$"),
        \(l) map(l, \(ln) suppressWarnings(loo::psis(rowSums(ln)))), 
        .names = "psis_{.col}"
      ), 
      across(
        starts_with("psis"), 
        lst(k = \(po) map_dbl(po, loo::pareto_k_values), lwt = \(po) map(po, \(pon) weights(pon, normalize = TRUE)[, 1])), 
        .names = "{.fn}_{.col}"
      ),
      across(matches("^(psis|lwt|k)"), lag),
    ) |> 
    rename_with(\(n) str_replace_all(
      n, 
      c(r"{log_lik_log_ratio}" = "log_ratio",
        r"{(k|lwt)_psis_patient(_.+)?_log_ratio}" = r"{\1\2}", 
        r"{^psis_patient(_.+)?_log_ratio}" = r"{psis\1}")
    )) |>   
    mutate(
      dplyover::across2(matches("^patient(_.+)?_log_lik$"), matches("^lwt(_.+)?"), psis_resample, .names = "approx_mean_{xcol}"),
      dplyover::across2(matches("^patient(_.+)?_log_lik_w$"), matches("^lwt(_+)?"), psis_resample, .names = "approx_mean_{xcol}"),
      across(matches("^(approx_)?mean"), \(m) map_dbl(m, \(mn) if (!is_null(mn)) sum(mn) else NA_real_), .names = "E_{.col}")
    ) |> 
    rename_with(\(n) str_replace(n, r"{E_(approx_)?mean}", r"{\1E}"))
}

redo_lfo_results <- function(lfo_res, lean = FALSE) {
  new_res <- lfo_res |> 
    lfo_drop_bad_approx() |> 
    arrange(n) |> 
    group_by(refit_n) |> 
    reframe(lfo_log_lik(first(fit), max_n = n())) |>   
    mutate(n = n + refit_n - 1)
  
  if (lean) {
    new_res <- new_res |>
      select(n, contains("E_"), k)
  }
  
  return(new_res)
}

clean_lfo_results <- function(lfo_res) {
  lfo_res |> 
    lfo_drop_bad_approx() |> 
    mutate(
      E_log_lik = if_else(is.na(k), E_patient_log_lik, approx_E_patient_log_lik),
      E_pfs_log_lik = if_else(is.na(k), E_patient_pfs_log_lik, approx_E_patient_pfs_log_lik),
      E_crcr_log_lik = if_else(is.na(k), E_patient_crcr_log_lik, approx_E_patient_crcr_log_lik),
      E_log_lik_w = if_else(is.na(k), E_patient_log_lik_w, approx_E_patient_log_lik_w),
      E_pfs_log_lik_w = if_else(is.na(k), E_patient_pfs_log_lik_w, approx_E_patient_pfs_log_lik_w),
      E_crcr_log_lik_w = if_else(is.na(k), E_patient_crcr_log_lik_w, approx_E_patient_crcr_log_lik_w),
    )
}

#' Bootstrap Expected Log Pointwise Predictive Density (ELPD) for Leave-Future-Out Cross-Validation
#'
#' This function performs bootstrapping to estimate the uncertainty in the Expected Log Pointwise 
#' Predictive Density (ELPD) for Leave-Future-Out (LFO) cross-validation results.
#'
#' @param lfo_res A data frame containing the results of LFO cross-validation, typically output 
#'        from the `lfo` function. Expected to contain columns with patient log-likelihoods and weights.
#' @param n_bootstrap Integer. The number of bootstrap samples to generate (default: 1000).
#'
#' @return A data frame with `n_bootstrap` rows, each containing a bootstrapped estimate of the ELPD.
#'         The columns correspond to different components of the log-likelihood (e.g., overall, PFS, CRCR).
#'
#' @details
#' The function performs the following steps:
#' 1. Removes any bad approximations from the LFO results.
#' 2. For each bootstrap iteration:
#'    a. Samples patients with replacement.
#'    b. Recalculates log-likelihoods and weights for the sampled patients.
#'    c. Computes the ELPD using PSIS (Pareto Smoothed Importance Sampling).
#' 3. Returns a data frame of bootstrapped ELPD estimates.
#'
#' This bootstrapping approach helps quantify the uncertainty in the ELPD estimate, 
#' which is crucial for model comparison and assessment in a time-series context, 
#' particularly for clinical trial data with multiple outcomes.
#'
lfo_bootstrap_elpd <- function(lfo_res, n_bootstrap = 1000) {
  get_patient_subset_col <- function(ln, i) {
    n_early_patients <- length(i) - ncol(ln)
    ln[, discard(i, \(x) x < n_early_patients) - n_early_patients]
  } 
  
  origin_res <- lfo_res |> 
    lfo_drop_bad_approx() |> 
    arrange(n)
  
  n_lfo_patients <- ncol(first(origin_res$patient_log_lik))
  
  map_dfr(seq(n_bootstrap), function(b) {
    bootstrap_i <- sample(n_lfo_patients, n_lfo_patients, replace = TRUE)
   
    origin_res |>
      mutate(
        across(matches("^patient(_pfs|_crcr)?_log_lik(_w)?$"), \(l)  map(l, \(ln) get_patient_subset_col(ln, bootstrap_i))),
        dplyover::across2(
          matches("^patient(_pfs|_crcr)?_log_lik$"), matches("^lwt(_pfs|crcr)?"), \(l, lw) psis_resample(l, lw, recalc_full = TRUE), .names = "mean_{xcol}"
        ),
        dplyover::across2(
          matches("^patient(_pfs|_crcr)?_log_lik_w$"), matches("^lwt(_pfs|crcr)?"), \(l, lw) psis_resample(l, lw, recalc_full = TRUE), .names = "mean_{xcol}"
        ),
      ) |> 
      transmute(across(matches("^mean"), \(m) map_dbl(m, \(mn) if (!is_null(mn)) sum(mn) else NA_real_), .names = "E_{.col}")) |> 
      rename_with(\(n) str_replace(n, r"{E_(approx_)?mean_patient}", r"{\1E}")) |> 
      summarize(across(everything(), sum))
  }) 
}

lfo_stacking_weights <- function(model_log_lik, log_lik_var = E_log_lik) {
  model_log_lik |> 
    map_dfr(clean_lfo_results, .id = "model") |> 
    select(model, n, {{ log_lik_var }}) |> 
    pivot_wider(names_from = model, values_from = {{ log_lik_var }}) |> 
    select(!n) |> 
    as.matrix() |> 
    loo::stacking_weights() |> 
    c() |> 
    set_names(names(model_log_lik))
} 

# nolint end: object_usage_linter