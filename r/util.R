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

plot_surv_ppc <- function(ppc_data, surv_interval_col, ic_col, rc_col, rep_surv_interval_col, exit_color = NULL) {
  ppc_data |> 
    mutate(
      max_week = max({{ surv_interval_col }} + {{ ic_col }}),
      usubjid = fct_reorder(usubjid, {{ surv_interval_col }} + {{ rc_col }} * max_week)
    ) |> 
    ggplot(aes(y = usubjid)) +
    stat_interval(aes(xdist = {{ rep_surv_interval_col }}), alpha = 0.5, linewidth = 1, .width = c(0.5, 0.8)) +
    geom_segment(
      aes(yend = usubjid, x = {{ surv_interval_col }}, xend = {{ surv_interval_col }} + {{ ic_col }}), linewidth = 0.25,
      data = \(d) filter(d, !{{ rc_col }})
    ) +
    geom_segment(
      aes(yend = usubjid, x = {{ surv_interval_col }}, xend = max_week), linewidth = 0.25, linetype = "dashed",
      data = \(d) filter(d, {{ rc_col }})
    ) +
    geom_point(aes(x = {{ surv_interval_col }}), size = 0.5) +
    geom_point(aes(x = {{ surv_interval_col }} + {{ ic_col }}, color = {{ exit_color }}), size = 0.5, data = \(d) filter(d, !{{ rc_col }})) +
    scale_color_discrete("", type = AZ_palette, label = c("FALSE" = "Non-response", "TRUE" = "Response"), aesthetic = c("color", "fill")) +
    # scale_color_ramp_discrete() +
    labs(y = "") +
    facet_grid(
      vars(trial), scales = "free_y", space = "free_y", switch = "y", labeller = labeller(trial = c("endometrial" = "Endometrial", "lung" = "Lung"))
    ) +
    theme(
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(), 
      panel.grid.minor.y = element_blank(), 
      legend.position = "top",
      strip.text.y.left = element_text(angle = 0)
    ) +
    guides(
      color_ramp = "none",  # Remove the interval legend
      color = guide_legend("")  # Keep only the confirmed_response legend
    ) +
    NULL
}