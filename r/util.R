#' Convert Kaplan-Meier estimates to a tibble (data frame) format 
#'
#' @param trt_data Analysis data 
#' @param key Identifier for the data group (e.g., treatment arm)
#' @param pfs_var Name of variable were PFS is stored in the data 
#'
#' @return tibble object with Kaplan-Meier results.
km_to_tibble <- function(trt_data, key, pfs_var, pfs_functions) { 
  with(
    base_prepare_pfs_stan_data(trt_data, pfs_var = pfs_var, pfs_functions), {
      interval_censored <- pfs_functions$identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure)[[1]]
      
      map_dfr(list(lb = pfs, ub = pfs + interval_censored), function(s) {
        pfs_functions$estimate_kaplan_meier(s, right_censored, max(s)) |>
          set_names(c("s", "n", "c", "e")) |>
          as_tibble() |> 
          mutate(t = seq(0, n() - 1))
      }, .id = "btype")
    }) |> 
    bind_cols(key)
}

get_km_res <- function(analysis_data, pfs_var, pfs_functions, ...) {
  analysis_data |>
    group_by(trial, ...) |>  
    group_map(\(trt_data, key) km_to_tibble(trt_data, key, pfs_var, pfs_functions), .keep = TRUE) |>  
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
    list_assign(!!!priors) %>% 
    list_assign(
      crcr_covar_effect_sd = rep(.$crcr_covar_effect_sd, .$n_covar),
      crcr_tumor_stim_pop_coef_sd = .$crcr_tumor_stim_pop_coef_sd[1:2],
    )
}

add_pfs_crcr_priors <- function(stan_data, crcr_priors, tumor_priors, pfs_priors) {
  add_confirmed_resp_priors(stan_data, crcr_priors) |> 
    list_assign(!!!tumor_priors, !!!pfs_priors) %>% 
    list_assign(
      covar_effect_sd = rep(.$covar_effect_sd, .$n_covar),
      tumor_stim_pop_coef_sd = .$tumor_stim_pop_coef_sd[1:2],
    )
}

create_pfs_crcr_initializer <- function(stan_data, n_causes = 2) {
  crcr_init_fun <- create_crcr_initializer(stan_data, n_causes)
  
  function(chain_id) {
    init_vals <- crcr_init_fun(chain_id) |> 
      list_assign(
        log_lambda_gp_intercept = with(stan_data, rnorm(1, log_lambda_gp_intercept_mean, log_lambda_gp_intercept_sd))
      )
    
    if (stan_data$add_trial_level) {
       init_vals <- init_vals |>  
        list_assign(
          log_lambda_gp_trial_intercept_sd = with(stan_data, abs(rnorm(1, sd = log_lambda_gp_trial_intercept_sd_sd)))
        )  
    }
    
    return(init_vals)
  }
}

# This function is used to generate a histogram of time-to-events for a single draw
sample_hist <- function(pred, breaks, ...) {
  # hist() is a base R function to generate histograms from data and provided breaks.
  hist(pmax(pmin(pred, max(breaks)), min(breaks)), breaks = breaks, plot = FALSE, ...)$count
}

# This function is used to treated_pfs_analysis_dataallow us to generate a distribution of histograms
rvar_sample_hist <- posterior::rfun(sample_hist, rvar_dots = FALSE)

