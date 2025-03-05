#' Sample from a Stan model and optionally save the output
#'
#' @param model A Stan model object
#' @param ... Additional arguments passed to the sample method
#' @param output_dir Directory to save output files
#' @param output_basename Base name for output files
#' @param timestamp Boolean, whether to include a timestamp in file names
#' @param no_save Boolean, if TRUE, don't save any output files
#'
#' @return A fitted Stan model object
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

add_confirmed_resp_priors <- function(stan_data, priors) {
  stan_data |> 
    list_assign(!!!priors) 
}

add_pfs_crcr_priors <- function(stan_data, crcr_priors, tumor_priors, pfs_priors) {
  add_confirmed_resp_priors(stan_data, crcr_priors) |> 
    list_assign(!!!tumor_priors, !!!pfs_priors)  
}

#' Generate a histogram of time-to-events for a single draw
#'
#' This function creates a histogram of time-to-event data for a single draw from a
#' posterior distribution. It uses R's base hist() function but returns only the counts,
#' not the full histogram object.
#'
#' @param pred A numeric vector of predicted time-to-event values
#' @param breaks A numeric vector specifying the breakpoints between histogram cells
#' @param ... Additional arguments passed to hist()
#'
#' @return A numeric vector of counts for each histogram bin
#'
sample_hist <- function(pred, breaks, ...) {
  # hist() is a base R function to generate histograms from data and provided breaks.
  hist(pmax(pmin(pred, max(breaks)), min(breaks)), breaks = breaks, plot = FALSE, ...)$count
}

# This function is used to treated_pfs_analysis_dataallow us to generate a distribution of histograms
rvar_sample_hist <- posterior::rfun(sample_hist, rvar_dots = FALSE)


Claude 3.5 Sonnet

2:33:38 pm

Certainly! Here's the documentation for the name_coef_indices function:

#' Name coefficient indices with meaningful labels
#'
#' This function takes a data frame with coefficient indices and adds meaningful labels
#' to these coefficients based on their index and the provided Stan data. It also
#' optionally adds trial labels.
#'
#' @param data A data frame containing the coefficient indices to be named
#' @param coef_idx_col The name of the column in 'data' that contains the coefficient indices
#' @param trial_col The name of the column in 'data' that contains trial identifiers (optional)
#' @param stan_data A list containing Stan data, including 'covar_design_matrix' and 'patient_trial'
#'
#' @return A modified data frame with additional columns:
#'   - 'covar': A factor column with meaningful names for each coefficient
#'   - 'trial': A factor column with trial labels (if trial_col is provided)
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

weeks_to_months <- function(weeks) weeks * 7 * 12 / 365.25
label_weeks_to_months <- scales::label_number(scale = weeks_to_months(1))
months_to_weeks <- function(months) months / weeks_to_months(1) 