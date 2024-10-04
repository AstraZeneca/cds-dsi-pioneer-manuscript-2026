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
