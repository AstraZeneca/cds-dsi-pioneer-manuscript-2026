sample_and_save <- function(model, ..., output_dir, output_basename, timestamp = TRUE, no_save = FALSE) {
  if (!no_save && !timestamp) {
    fit <- model$sample(..., output_dir = output_dir, output_basename = output_basename)
  } else { 
    fit <- model$sample(...)
    
    if (!no_save) {
      fit$save_output_files(dir = output_dir, basename = output_basename, random = FALSE, timestamp = timestamp)
      fit$save_profile_files(dir = output_dir, basename = output_basename, random = FALSE, timestamp = timestamp)
    }
  }
  
  return(fit)
}

#' Convert Kaplan-Meier estimates to a tibble (data frame) format 
#'
#' @param trt_data Analysis data 
#' @param key Identifier for the data group (e.g., treatment arm)
#' @param pfs_var Name of variable were PFS is stored in the data 
#'
#' @return tibble object with Kaplan-Meier results.
km_to_tibble <- function(trt_data, key, pfs_var) { 
  stan_data <- base_prepare_pfs_stan_data(trt_data, pfs_var = pfs_var) |> 
    magrittr::extract(c("pfs", "interval_censored", "right_censored"))
  
  lst(
    lb = survfit2(Surv(pfs + 1, 1 - right_censored) ~ 1, stan_data),
    ub = survfit2(Surv(pfs + interval_censored + 1, 1 - right_censored) ~ 1, stan_data),
  ) |> 
    map_dfr(broom::tidy, .id = "btype") |>  
    select(t = time, s = estimate, n = n.risk, c = n.censor, e = n.event, btype) |> 
    bind_cols(key)
}

get_km_res <- function(analysis_data, pfs_var, ...) {
  analysis_data |>
    group_by(trial, ...) |>  
    group_map(\(trt_data, key) km_to_tibble(trt_data, key, pfs_var), .keep = TRUE) |>  
    bind_rows() 
} 

cmdstan_expose_pfs_functions <- function(util_file, pfs_functions_file) {
  pseudo_model_code <- paste(c("functions {", read_file(util_file), read_file(pfs_functions_file), "}"), collapse="\n")
  functions_hash <- rlang::hash(pseudo_model_code)
  model_name <- paste0("pfs-functions-", functions_hash)
  ## note: cmdstanr somehow only compiles standalone functions
  ## whenever one is compiling the model (and not allowing to export
  ## the functions if one is not compiling it). This is why
  ## force_compile=TRUE is a save option
  ##pseudo_model <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(pseudo_model_code), compile_standalone=TRUE, force_compile=TRUE, stanc_options=list(name=paste0("model-functions-", functions_hash)))
  ##pseudo_model$functions
  ## but things seem to work ok if we abuse a bit the internals... tested with cmdstanr 0.6.1
  ## note that we have to set the model name manually to a
  ## determinstic string (depending only on the stan functions being
  ## compiled)
  stan_file <- cmdstanr::write_stan_file(pseudo_model_code)
  pseudo_model <- cmdstanr::cmdstan_model(stan_file, stanc_options=list(name=model_name))
  pseudo_model$functions$existing_exe <- FALSE
  pseudo_model$functions$external <- FALSE
  stancflags_standalone <- c("--standalone-functions", paste0("--name=", model_name))
  pseudo_model$functions$hpp_code <- cmdstanr:::get_standalone_hpp(stan_file, stancflags_standalone)
  pseudo_model$expose_functions(FALSE, FALSE) ## will return the functions in an environment
  pseudo_model$functions
}

add_confirmed_resp_priors <- function(stan_data, priors) {
  stan_data |> 
    list_assign(!!!priors) 
    # list_assign(
    #   crcr_tumor_stim_pop_coef_sd = .$crcr_tumor_stim_pop_coef_sd,
    # )
}

add_pfs_crcr_priors <- function(stan_data, crcr_priors, tumor_priors, pfs_priors) {
  add_confirmed_resp_priors(stan_data, crcr_priors) |> 
    list_assign(!!!tumor_priors, !!!pfs_priors)  
    # list_assign(
    #   tumor_stim_pop_coef_sd = .$tumor_stim_pop_coef_sd[1:2],
    # )
}

# This function is used to generate a histogram of time-to-events for a single draw
sample_hist <- function(pred, breaks, ...) {
  # hist() is a base R function to generate histograms from data and provided breaks.
  hist(pmax(pmin(pred, max(breaks)), min(breaks)), breaks = breaks, plot = FALSE, ...)$count
}

# This function is used to treated_pfs_analysis_dataallow us to generate a distribution of histograms
rvar_sample_hist <- posterior::rfun(sample_hist, rvar_dots = FALSE)

name_coef_indices <- function(data, coef_idx_col, trial_col, stan_data) {
  data |> 
    mutate(
      covar = case_when(
        {{ coef_idx_col }} == 1 ~ "baseline sum of tumor sizes",
        {{ coef_idx_col }} == 2 ~ "first post-treatment sum of tumor sizes",
        {{ coef_idx_col }} - 2 <= ncol(stan_data$covar_design_matrix) ~ 
          colnames(stan_data$covar_design_matrix)[pmax(1, {{ coef_idx_col }} - 2)] |> 
          str_replace(r"{factor\((.+),\sordered\s=\sFALSE\)}", "\\1 "),
        TRUE ~ "confirmed response"
      ) |> as_factor(),
      trial = if(!is_null(trial_col)) factor({{ trial_col }}, labels = levels(stan_data$patient_trial)),
    )
}

get_fake_stan_data_list <- function(prior_res, origin_stan_data, n = 5) {
  get_all_confirmed_response(prior_res) |> 
    select(starts_with("rep_")) |> 
    unnest_rvars() |> 
    filter(.draw <= n) |> 
    rename_with(\(n) str_remove(n, "^rep_")) |>
    select(!c(.chain, .iteration)) |> 
    group_by(.draw) |>
    group_map(\(d, k, ...) list_assign(origin_stan_data, !!!d, draw = first(k$.draw), confirmed_response_interval_censored = rep(0, nrow(d)))) 
}

# recist_response <- function(baseline_sum, current_sum) {
#   if (!is.numeric(baseline_sum) || !is.numeric(current_sum) || 
#       baseline_sum <= 0 || current_sum < 0) {
#     stop("Inputs must be positive numbers, with baseline > 0")
#   }
#   
#   absolute_change <- current_sum - baseline_sum
#   percent_change <- absolute_change / baseline_sum 
#   
#   case_when(
#     current_sum == 0 ~ "CR",
#     percent_change <= -0.3 ~ "PR",
#     percent_change >= 0.2 & absolute_change >= 5 ~ "PD",
#     TRUE ~ "SD"
#   )
# }

#' Determine RECIST 1.1 Response
#'
#' This function calculates the RECIST 1.1 response category based on measurements
#' of target lesions, and optionally considers non-target lesions and new lesions.
#'
#' @param baseline_sld Baseline sum of longest diameter
#' @param current_sld Current sum of longest diameter
#' @param nadir_sld The smallest sum of measurements observed so far (default is NULL, which means baseline is used as nadir)
#' @param include_non_target Logical indicating whether to consider non-target lesions (default is FALSE)
#' @param non_target_response A character string: "CR", "SD", or "PD" (only used if include_non_target = TRUE)
#' @param new_lesions Logical indicating whether new lesions have appeared (default is FALSE)
#'
#' @return A character string indicating the RECIST 1.1 response category
#'
determine_recist_response <- function(baseline_sld, current_sld, nadir_sld = NULL, include_non_target = FALSE, non_target_response = "NON-CR/NON-PD", new_lesions = FALSE) {
  # Input validation
  if (!is.numeric(baseline_sld) || !is.numeric(current_sld) || baseline_sld < 0 || current_sld < 0) {
    stop("baseline_sld and current_sld must be non-negative numeric values")
  }
  if (!is.null(nadir_sld) && (!is.numeric(nadir_sld) || nadir_sld < 0)) {
    stop("nadir_sld must be a non-negative numeric value or NULL")
  }
  if (!is.logical(include_non_target) || !is.logical(new_lesions)) {
    stop("include_non_target and new_lesions must be logical values")
  }
  if (!is.na(non_target_response) && !non_target_response %in% c("CR", "NON-CR/NON-PD", "PD", "NE")) {
    stop("non_target_response must be 'CR', 'NON-CR/NON-PD', 'PD', or NA, got: ", non_target_response)
  }
  
  # If nadir sum not provided, use baseline as nadir
  if (is_null(nadir_sld)) {
    nadir_sld <- min(baseline_sld, current_sld)
  } else {
    nadir_sld <- min(nadir_sld, current_sld)  # Update nadir if current sum is smaller
  }
  
  # Calculate changes
  change_from_baseline <- (current_sld - baseline_sld) / baseline_sld
  absolute_diff_from_nadir <- current_sld - nadir_sld
  change_from_nadir <- absolute_diff_from_nadir / nadir_sld
  
  # Define progression for target lesions (≥20% increase from nadir AND ≥5mm absolute increase)
  target_progression <- current_sld > nadir_sld && change_from_nadir >= 0.2 && absolute_diff_from_nadir >= 5 
  
  # If considering only target lesions
  if (!include_non_target && !new_lesions) {
    case_when(
      current_sld == 0 ~ "CR",
      change_from_baseline <= -0.3 ~ "PR",
      target_progression ~ "PD",
      .default = "SD"
    )
  } else {
    case_when(
      # If including non-target lesions or new lesions
      target_progression || (!is.na(non_target_response) && non_target_response == "PD") || new_lesions ~ "PD",
      current_sld == 0 && (non_target_response == "CR" || is.na(non_target_response)) ~ "CR",
      change_from_baseline <= -0.3 ~ "PR",
      .default = "SD"
    )
  }
}

weeks_to_months <- function(weeks) weeks * 7 * 12 / 365.25
label_weeks_to_months <- scales::label_number(scale = weeks_to_months(1))
months_to_weeks <- function(months) months / weeks_to_months(1) 

lognormal_sd <- function(mu = 0, sigma) {
  # Calculate standard deviation of lognormal variable X
  # where log(X) ~ N(mu, sigma)
  #
  # Args:
  #   mu: mean parameter of the normal distribution in log space
  #   sigma: standard deviation parameter of the normal distribution in log space
  #
  # Returns:
  #   standard deviation of the lognormal random variable X
  
  sqrt((exp(sigma^2) - 1) * exp(2*mu + sigma^2))
} 
