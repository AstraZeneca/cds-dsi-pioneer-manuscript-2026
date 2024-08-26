#' Convert Kaplan-Meier estimates to a tibble (data frame) format 
#'
#' @param trt_data Analysis data 
#' @param key Identifier for the data group (e.g., treatment arm)
#' @param pfs_var Name of variable were PFS is stored in the data 
#'
#' @return tibble object with Kaplan-Meier results.
km_to_tibble <- function(trt_data, key, pfs_var, tumor_priors, pfs_priors, pfs_functions) { 
  with(
    prepare_pfs_stan_data(trt_data, tumor_priors, pfs_priors, pfs_var = pfs_var, pfs_functions), {
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

get_km_res <- function(analysis_data, pfs_var, tumor_priors, pfs_priors, pfs_functions) {
  analysis_data |>
    group_by(trial, treated) |>  
    group_map(\(trt_data, key) km_to_tibble(trt_data, key, pfs_var, tumor_priors, pfs_priors, pfs_functions), .keep = TRUE) |>  
    bind_rows() 
} 

#' Function to reorganize t tumor size measure arrays.
#' 
#' Information about the number of measures and the intervals of measurement of tumor sizes are typically 
#' stored in a flat array due to the lack of support for ragged arrays in Stan. This function converts this information
#' into a structure easier to use in R. 
#'
#' @param measures The intervals in which each tumor measure is done. 
#' @param n_measures The number of measures done for each tumor.
#' @param n_patient_tumors The number of tumors per patient. 
#'
#' @return Nest list
list_measures <- function(measures, n_measures, n_patient_tumors) { 
  split(measures, rep(seq_along(n_measures), n_measures)) |>  
    split(rep(seq_along(n_patient_tumors), n_patient_tumors))
}

#' Drop particular measures from the the tumor size data.
#' 
#' This is used to generate simulated data with interval censoring. 
#'
#' @param settings Stan configurations. 
#' @param measures Measures to drop. 
#' @param keep_only Instead of dropping measures, keep only the ones specified.
#'
#' @return Updated Stan data with specified measures dropped.
drop_missing_measures <- function(settings, measures = NULL, keep_only = FALSE) {
  if (!is_null(measures) && length(measures) > 0) {
    updated_settings <- if (is.list(measures)) {
      if (keep_only) {
        settings %>%  
          list_modify(
            t_measure = measures,
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map2(measures, \(t, m) keep_at(t, m - min(m) + 1)) |>  
              unlist()
          ) 
      } else {
        settings %>%  
          list_modify(
            t_measure = list_measures(.$t_measure, .$n_measures - 1) |> 
              map2(measures, \(t, m) setdiff(t, m)),
            tumor_size = list_measures(.$tumor_size, rep(.$n_measures, .$n_patient_tumors)) |> 
              map2(rep(measures, .$n_patient_tumors), \(t, m) discard_at(t, m + 1)) |>  
              unlist()
          ) 
      }
    } else {
      if (keep_only) {
        settings %>%  
          list_modify(
            t_measure = map(seq(sum(.$n_patient_tumors)), \(i) measures),
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) keep_at(t, measures - min(measures) + 1)) |>  
              unlist()
          ) 
      } else {
        settings %>%  
          list_modify(
            t_measure = with(., list_measures(t_measure, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) setdiff(t, measures)),
            tumor_size = with(., list_measures(tumor_size, n_measures, n_patient_tumors)) |> 
              list_flatten() |> 
              map(\(t) discard_at(t, measures - min(measures) + 1)) |>  
              unlist()
          ) 
      }
    }
    
    updated_settings %>%
      list_assign(
        n_measures = map_int(.$t_measure, length),
        t_measure = unlist(.$t_measure),
      ) 
  } else {
    settings
  }
}

#' ORR Model Stan initializer factory. 
#'
#' @param stan_data Analysis data in list form for Stan. 
#'
#' @return Initialization function.
create_orr_initializer <- function(stan_data) {
  function(chain_id) { 
    init_vals <- lst(
    )
    
    if (stan_data$add_trial_level) {
      init_vals <- init_vals |> 
        list_assign(
          raw_log_lambda_gp_trial_intercept = if (stan_data$add_trial_level) rnorm(stan_data$n_trials),
          log_lambda_gp_trial_intercept_sd = abs(rnorm(1, 0, stan_data$log_lambda_gp_trial_intercept_sd_sd)),
        )
    }
    
    if (stan_data$add_tumor_location_level) {
      init_vals <- init_vals |> 
        list_assign(
        )
    }
    
    return(init_vals)
  }
}

#' Tumor size model Stan initializer factory. 
#'
#' @param stan_data Analysis data in list form for Stan. 
#'
#' @return Initialization function.
create_pfs_initializer <- function(stan_data) {
  function(chain_id) { 
    n_covar <- with(stan_data, if_else(tumor_hazard_type > 0 && tumor_hazard_type != 4, 
                                       if_else(tumor_hazard_type < 3, 
                                               tumor_hazard_type + 1, 
                                               if_else(tumor_hazard_type == 5, 5, 1)), 
                                       0))
    
    init_vals <- lst(
      tumor_stim_pop_intercept = abs(rnorm(1, 0, stan_data$tumor_stim_pop_intercept_sd)),
      tumor_stim_pop_coef = abs(rnorm(n_covar, 0, stan_data$tumor_stim_pop_coef_sd[1:n_covar])), 
    )
    
    if (stan_data$add_trial_level) {
      init_vals <- init_vals |> 
        list_assign(
          raw_log_lambda_gp_trial_intercept = if (stan_data$add_trial_level) rnorm(stan_data$n_trials),
          raw_tumor_stim_trial_coef = map(seq(stan_data$n_trials), \(...) rnorm(n_covar + 1)),
          log_lambda_gp_trial_intercept_sd = abs(rnorm(1, 0, stan_data$log_lambda_gp_trial_intercept_sd_sd)),
          tumor_stim_trial_coef_sd = abs(rnorm(n_covar + 1, 0, stan_data$tumor_stim_trial_coef_sd_sd)),
        )
    }
    
    if (stan_data$add_tumor_location_level) {
      init_vals <- init_vals |> 
        list_assign(
          raw_tumor_stim_location_coef = map(seq(stan_data$n_tumor_locations), \(...) rnorm(n_covar + 1)),
          tumor_stim_location_coef_sd = abs(rnorm(n_covar + 1, 0, stan_data$tumor_stim_location_coef_sd_sd)),
        )
    }
    
    return(init_vals)
  }
}

#' Run Stan sampling on given simulation data. 
#'
#' @param settings Stan data/settings 
#' @param d Analysis data 
#' @param max_measures deprecated setting 
#' @param ... Any other parameters to pass to cmdstanr::sample().
#' @param gen_pfs Should the model generated simulated data. 
#' @param ignore_interval_censoring The model should ignore interval censoring. 
#' @param drop_measures deprecated setting 
#' @param keep_only deprecated setting 
#' @param no_init deprecated setting
#'
#' @return cmdstanr fit object
fit_sim_data <- function(
  settings, d, max_measures, ..., gen_pfs = TRUE, ignore_interval_censoring = FALSE, drop_measures = NULL, keep_only = FALSE, no_init = FALSE
) { 
  settings <- settings |> 
    list_modify(
      gen_pfs = gen_pfs,
      fit_data = TRUE,
      pfs = d$pfs, 
      right_censored = d$right_censored,
      ignore_interval_censoring = ignore_interval_censoring
    )
  
  settings |> 
    pfs_model$sample(
      refresh = 0,
      parallel_chains = 4,
      init = create_pfs_initializer(settings),
      ...
    )
}

#' Identify the patient and tumor for each posterior tumor size time-series. 
#'
#' @param rv rvar data for tumor sizes 
#' @param settings Stan data
#'
#' @return rvar data set but with patient_id, tumor_index, and t
prepare_tumor_size_rvars <- function(rv, patient_measures, t_offset, settings) {
  rv |> 
    ungroup() |> 
    mutate( # A bunch of acrobatics to get the IDs, indices, and intervals right.
      patient_id = with(settings, rep(1:n_patients, n_patient_tumors * patient_measures)),
      tumor_index = ((tumor_index - 1) %/% patient_measures) + 1,
      t = rep(0:(patient_measures - 1), sum(settings$n_patient_tumors)) - t_offset
    ) |>
    mutate(tumor_index = tumor_index - min(tumor_index) + 1, .by = patient_id) |>
    rename(tumor_size = rep_tumor_size)
}

spread_param_rvars <- function(f, ...) { 
  f |> 
    spread_rvars(
      tumor_stim_pop_intercept, tumor_stim_pop_coef[t],
      log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho,
      base_cond_expected_pfs, one_tumor_cond_expected_pfs, base_cond_median_pfs, one_tumor_cond_median_pfs,
      ...
    ) |> 
    ungroup() |> 
    pivot_wider(names_from = t, values_from = tumor_stim_pop_coef, names_prefix = "tumor_stim_pop_coef_")
}

get_pfs_sim_seed_draws <- function(stan_data, num_sim, pfs_model) {
  pfs_res <- pfs_model$sample(data = stan_data, refresh = 0, parallel_chains = 4)
  
  pfs_res |> 
    spread_draws(rep_pfs[patient_index], rep_right_censored[patient_index]) |> 
    ungroup() |>
    # The PFS from each draw will be used as a simulation dataset 
    filter(.draw <= num_sim) |> 
    select(.draw, pfs = rep_pfs, right_censored = rep_right_censored) |> 
    nest(sim_data = !.draw) |>
    left_join( # Get the parameters that generated that data
      pfs_res |> 
        spread_param_rvars() |> 
        unnest_rvars() |> 
        select(!c(.iteration, .chain)) |> 
        pack(true = !.draw),
      by = ".draw"
    )
}

run_simulation <- function(
    sim_id, stan_data, d, lstm, output_dir, output_name, ignore_interval_censoring, gen_pfs, pfs_model, 
    .thin = NULL, .ndraws = NULL, keep_fit = FALSE
) { 
  d <- mutate(d, map2_dfr(pfs, lstm, function(s, m) { 
      m_union <- reduce(m, \(a, n) union(a, n)) 
      list(pfs = m_union |> discard(\(t) t > s) |> max(), right_censored = s >= max(m_union))
  }))
  
  fit <- stan_data |> 
    list_modify(
      gen_pfs = gen_pfs,
      fit_data = TRUE,
      pfs = d$pfs, 
      right_censored = d$right_censored,
      ignore_interval_censoring = ignore_interval_censoring
    ) %>% 
    pfs_model$sample(
      refresh = 0,
      parallel_chains = 4,
      init = create_pfs_initializer(.),
      output_basename = str_glue("{output_name}_{sim_id}"),
      output_dir = output_dir,
      thin = .thin,
    )
 
  if (!is_null(.ndraws)) { 
    res <- fit |> # Get posterior draws from simulation fit. 
      spread_param_rvars(ndraws = .ndraws) |> 
      pack(est = everything())
  } else {
    res <- fit |> # Get posterior draws from simulation fit. 
      spread_param_rvars() |> 
      pack(est = everything())
  }
  
  if (keep_fit) {
    res <- res |> 
      mutate(sim_fit = list(fit))
  }
  
  return(res)
}

#' Run simulation fit.
#' 
#' Initially, used for simulation-based calibration but can also be used for any analysis of fake data. 
#'
#' @param pfs_model cmdstanr model object. 
#' @param stan_data Analysis and configuration data formatted for use in Stan model. 
#' @param num_sim Number of simulations to carry out.
#' @param output_name Stan samples output file name prefix 
#' @param output_dir Where to store Stan samples output files 
#' @param ignore_interval_censoring Do not adjust for interval censoring 
#' @param keep_fit Keep cmdstanr fit object in function output. 
#' @param gen_pfs Should the simulation fit runs also generate simulated PFS.
#' @param .thin Rate of thinning in output samples (needed for SBC).
#' @param .ndraws Number of draws to extract from posterior. 
#' @param reuse_data Do not generate new simulation.
#' @param drop_measures Which measure intervals to drop to simulate interval censoring. 
#' @param keep_only Keep only the measure intervals specified and drop all others.
#' @param ... Any other parameters to pass to fit_sim_data() .
#'
#' @return Data frame with all the simulation posterior samples (and possibly fit object).
run_sbc_sims <- function(
  pfs_model, stan_data, num_sim, output_name,  output_dir = file.path(tmp_dir, "fit"), 
  ignore_interval_censoring = FALSE,keep_fit = FALSE, gen_pfs = FALSE, .thin = 4, .ndraws = 1000,
  reuse_data = NULL, drop_measures = NULL, keep_only = FALSE, ... 
) {
  
  stan_data <- stan_data |> 
    drop_missing_measures(drop_measures, keep_only)
  
  lstm <- with(stan_data, list_measures(t_measure, n_measures, n_patient_tumors))
 
  all_sim_data <- if (is_null(reuse_data)) { 
    get_pfs_sim_seed_draws(stan_data, num_sim, pfs_model)
  } else {
    reuse_data |> select(.draw, sim_data, true)
  }
  
  all_sim_data |> 
    mutate(# For each simulation dataset fit the model and extract the posteriors for each of the model parameters 
      # sim_data = map(
      #   sim_data, 
      #   \(d) mutate(d, 
      #               map2_dfr(pfs, lstm, function(s, m) { 
      #                 m_union <- reduce(m, \(a, n) union(a, n)) 
      #                 list(pfs = m_union |> discard(\(t) t > s) |> max(), right_censored = s >= max(m_union))
      #               })
      #   )),
                              
      furrr::future_pmap_dfr(.progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
      # pmap_dfr(
        lst(sim_id = .draw, sim_data, lstm),
        rum_simulation,
        output_dir = output_dir, output_name = output_name, ignore_interval_censoring = ignore_interval_censoring, 
        gen_pfs = gen_pfs, .thin = .thin, .ndraws = .draws
      ),
    ) 
}

#' Plot Kaplan-Meier survival curves 
#'
#' @param rvars_data Posterior samples data 
#' @param analysis_data Analysis data used for fit 
#' @param facet_arm Facet by `arm_var`, otherwise by trial 
#' @param group_arm Group output by `group_var`, otherwise no grouping. 
#' @param arm_var Data column to facet by. 
#' @param group_var Data column to group by. 
#' @param model_labels Labels use for models in legend. 
#'
#' @return ggplot2 plot
plot_km_rvars <- function(rvars_data, analysis_data, facet_arm = FALSE, group_arm = FALSE,
                          arm_var = arm, group_var = arm,
                          model_labels = c("no_tumor" = "Baseline Model", "tumor_change" = "Proportional Change Model", "two_tumor" = "Linear Model")) {
  plot_obj <- rvars_data |>  
    ggplot() +
    labs(x = "t [Week]", y = latex2exp::TeX("$S(t)$")) +
    NULL
  
  if (!facet_arm) {
    plot_obj <- plot_obj + facet_wrap(vars(trial), labeller = labeller("{{arm_var}}" := model_labels)) 
  } else {
    plot_obj <- plot_obj + facet_grid(vars({{arm}}), vars(trial), labeller = labeller("{{arm_var}}" := model_labels)) 
  }
  
  if (group_arm) {
    plot_obj +
      stat_lineribbon(aes(x = t - 1, ydist = trial_km_est, fill = {{group_var}}, color = {{group_var}}), step = TRUE, alpha = 0.5, .width = 0.8) +
      scale_color_discrete("", type = AZ_palette, aesthetics = c("color", "fill"))
  } else {
    plot_obj +
      stat_lineribbon(aes(x = t - 1, ydist = trial_km_est), fill = AZ_platinum, step = TRUE, alpha = 0.25, .width = 0.8) 
  }
}

#' Plot Kaplan-Meier survival curves for multiple fits 
#'
#' @param data_list List of posterior samples' data 
#' @param analysis_data Analysis data used for fit 
#' @param ... Other parameters to pass to `plot_km_rvars` 
#'
#' @return ggplot2 plot object.
plot_km <- function(data_list, analysis_data, ...) {
  data_list |> 
    map(\(r) recover_types(r, select(analysis_data, trial))) |> 
    map_dfr(\(r) spread_rvars(r, trial_km_est[trial, t]), .id = "arm") |>
    plot_km_rvars(analysis_data, arm_var = arm, ...)
}


# This function is used to generate a histogram of time-to-events for a single draw
sample_hist <- function(pred, breaks, ...) {
  # hist() is a base R function to generate histograms from data and provided breaks.
  hist(pmax(pmin(pred, max(breaks)), min(breaks)), breaks = breaks, plot = FALSE, ...)$count
}

# This function is used to treated_pfs_analysis_dataallow us to generate a distribution of histograms
rvar_sample_hist <- posterior::rfun(sample_hist, rvar_dots = FALSE)

#' Produce a probabilistic histogram from rvar samples
#' 
#' The difference between a probabilistic histogram and a regular histogram is that it shows the uncertainty about the distribution. 
#'
#' @param data_list List of cmdstanr fits 
#' @param stan_data Data used for fitting model 
#' @param hist_breaks Histogram breaks 
#'
#' @return ggplot2 plot object
plot_pfs_hist_posterior <- function(data_list, stan_data, hist_breaks = seq(10, 150, 10)) {
  model_labels = c("no_tumor" = "Baseline Model", "tumor_change" = "Proportional Change Model", "two_tumor" = "Linear Model")
  
  plot_obj <- data_list |> 
    map_dfr(\(f) spread_rvars(f, rep_pfs[i], rep_right_censored[i]) |> mutate(trial = stan_data$patient_trial), .id = "arm") |> 
    group_by(arm, trial) |> 
    # bindist is the distribution of histogram size at each bin.
    reframe(t = hist_breaks, bindist = rvar_sample_hist(rep_pfs, c(0, hist_breaks))) |> 
    filter(t < max(t)) %>% 
    bind_rows(
      group_by(., arm, trial) %>%
        filter(t %in% range(t)) |> 
        mutate(t = c(0, max(t) + min(hist_breaks)))
    ) |> 
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = bindist),
                    color = AZ_platinum, fill = AZ_platinum, alpha = 0.25, linewidth = 2,
                    step = "mid",
                    .width = c(0.5, 0.8), show.legend = FALSE) + 
    scale_x_continuous("t [Week]", breaks = seq(0, max(hist_breaks) + min(hist_breaks), 20)) +
    NULL
 
  if (length(data_list) > 1) {
    plot_obj + facet_grid(vars(arm), vars(trial), labeller = labeller(arm = model_labels)) 
  } else {
    plot_obj + facet_wrap(vars(trial), scales = "free_y") 
  }
}

#' Plot distribution of baseline hazard 
#'
#' @param data_list List of cmdstanr fit objects 
#' @param analysis_data Analysis data used 
#' @param ci_width Credible interval widths 
#'
#' @return ggplot2 plot object
plot_baseline_hazard <- function(data_list, analysis_data, ci_width = 0.8) {
  data_list |> 
    map(\(f) recover_types(f, select(analysis_data, trial))) |> 
    map_dfr(\(f) spread_rvars(f, log_trial_lambda[trial, t]), .id = "fit_type") |> 
    mutate(trial_lambda = exp(log_trial_lambda)) |> 
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = trial_lambda, color = fit_type, fill = fit_type, alpha = fit_type), step = TRUE, .width = ci_width) +
    geom_rug(aes(week), alpha = 0.125, 
             data = analysis_data |> 
               transmute(
                 trtp, trial, 
                 patient_visits = map(patient_tumors, \(tu) unnest(tu, tumor_history) |> distinct(week))
               ) |> 
               unnest(patient_visits)) +
    scale_color_viridis_d("", aesthetics = c("fill", "color"), labels = c(prior = "Prior", treated = "Posterior")) +
    scale_alpha_manual("", values = c(prior = 0.125, treated = 0.5)) +
    facet_wrap(vars(trial)) +
    labs(title = "Baseline hazard", y = latex2exp::TeX(r"{$\lambda_{st}$}"), caption = "Ribbons shown are for the 80% CI.\nRugs below x-axis show the distribution of assessment weeks.") +
    guides(alpha = "none") +
    theme(legend.position = "bottom") +
    NULL
}

#' Estimate out-of-sample loss 
#'
#' @param fit cmdstanr fit object 
#' @param fit_loo leave-one-out object 
#' @param stan_data Data used in the model 
#' @param loss_fn Loss function 
#' @param insample Run in-sample instead 
#'
#' @return Pointwise out-of-sample loss
estimate_oos_loss <- function (fit, fit_loo, stan_data, loss_fn, insample = FALSE) {
  log_lik <- fit$draws("log_lik", format = "matrix")
  
  observed_pfs <- with(stan_data, map2(pfs, interval_censored, \(p, i) runif(nrow(log_lik), p, p + i))) 
  rep_draws <- map(c("rep_pfs", "rep_right_censored"), \(v) fit$draws(v, format = "matrix")) |> 
    map(\(r) split(r, rep(seq(ncol(r)), each = nrow(r)))) |> 
    set_names(c("pfs", "right_censored"))
  
  calculated_loss <- lst(observed_pfs, observed_right_censored = stan_data$right_censored, !!!rep_draws) |> 
    unname() |> 
    pmap(loss_fn) |> 
    simplify2array()
  
  if (insample) {
    calculated_loss |> 
      plyr::aaply(2, mean) |> 
      enframe(name = NULL) |> 
      mutate(trial = stan_data$patient_trial, id = seq(n()))
  } else {
    calculated_loss |> 
      E_loo(fit_loo$psis_object, type = "mean", log_ratios = -log_lik) |> 
      as_tibble() |> 
      mutate(trial = stan_data$patient_trial, id = seq(n()))
  }
}

#' Plot estimation loss 
#'
#' @param data calculated loss data 
#' @param binwidth Histogram bin width 
#'
#' @return ggplot2 plot object
plot_loss <- function(data, binwidth) {
  data |> 
    ggplot(aes(value)) +
    stat_histinterval(
      aes(fill = model, color = model), point_interval = mean_qi, .width = c(0.5, 0.8), alpha = 0.5, breaks = ggdist::breaks_fixed(width = binwidth)
    ) +
    labs(x = "Loss", y = "") +
    facet_grid(vars(trial), vars(loss_type), scales = "free", margins = "trial") +
    theme(legend.position = "top", strip.text.y = element_text(angle = 0)) +
    NULL
}

#' Provide credible interval summary for hazard ratios are a specified grid representing tumor size changes 
#'
#' @param fit cmdstanr fit  
#' @param var Hazard ratio parameter name in model 
#' @param tumor_size_pairs Wide tumor size data 
#' @param stan_data Data used in model 
#' @param diff_var Calculate different in hazard ratio from the population hazard ratio 
#' @param ci_width The credible interval widths 
#' @param no_summary Just return the rvar data 
#'
#' @return Tumor-wise size time-series and time-invariant hazard ratio
get_tumor_hazard_ratio_summary <- function(fit, var, tumor_size_pairs, stan_data, diff_var = FALSE, ci_width = c(0, 0.5, 0.8), no_summary = FALSE) {
  grid_ids <- tibble(id = stan_data$grid_tumors) |> 
    mutate(n = seq(n()))
  
  tumor_size_pairs <- tumor_size_pairs |> 
    mutate(id = seq(n())) |> 
    semi_join(grid_ids, by = "id") |> 
    rename(tumor_trial = trial)
 
  # Translate CI widths to quantile() style probabilities 
  p_data <- tibble(
    p = (ci_width / 2) %>% { c(0.5 - ., 0.5 + .) } |> unique(),
  ) |> 
    mutate(pid = seq(n()))
  
  ratio_rvars <- fit |> 
    gather_rvars({{ var }}) |>
    left_join(grid_ids, by = "n", relationship = "many-to-one") |> 
    left_join(tumor_size_pairs, by = "id", relationship = "many-to-one") 
  
  if (diff_var) { # We want to see the difference between the hazard ratio from the baseline population level.
    ratio_rvars <- ratio_rvars |> 
      left_join(gather_rvars(fit, pop_tumor_effect[n]) |> select(base_value = .value, n), by = "n", relationship = "many-to-one") |> 
      mutate(.value = .value - base_value)
  }
  
  ratio_rvars <- ratio_rvars |> 
    mutate(.value = exp(.value)) 
  
  if (no_summary) {
    return(ratio_rvars)
  } else {
    ratio_rvars |> 
      rowwise() |> 
      # Turns out this is faster than median_qi()
      mutate(.value = list(enframe(quantile(.value, probs = p_data$p), name = "pid", value = "q"))) |> 
      ungroup() |>
      unnest(.value) |> 
      left_join(p_data, by = "pid") # I need to do this because quantile.rvar just adds indices not the p's
  }
}

get_sim_tumor_stan_data <- function(n_patients, patient_measures, t_offset, tumor_priors) {
  lst(
    fit_tumor_data = FALSE, # Prior prediction
    gen_tumor_sizes = TRUE,
    predict_missing_sizes = FALSE,
    multilevel_patient = FALSE,
    multilevel_tumor = FALSE,
    
    n_patients,
    n_trials = 4,
    patient_trial = sample(4, n_patients, replace = TRUE), 
    n_patient_tumors = 1 + rbinom(n_patients, 5, 0.4), # We're not yet generatively modeling the number of tumors
    n_measures = rep(patient_measures, sum(n_patient_tumors)),
    t_measure = rep(0:(patient_measures - 1), sum(n_patient_tumors)) - t_offset,
    tumor_size = rep(0, sum(n_measures)),
    tumor_location = sample(1:50, sum(n_patient_tumors), replace = TRUE), 
    n_tumor_locations = n_distinct(tumor_location),
    
    !!!tumor_priors  
  )
}

get_fake_tumor_data <- function(res, patient_measures, t_offset, stan_data) { 
  res |> 
    spread_rvars(rep_tumor_size[tumor_index], ndraws = 1) |> 
    prepare_tumor_size_rvars(patient_measures, t_offset, stan_data) |> 
    unnest_rvars()
}

get_tumor_ppc_draws <- function(res, patient_measures, t_offset, stan_data) {
  res |> 
    spread_rvars(rep_tumor_size[tumor_index]) |>
    ungroup() |> 
    prepare_tumor_size_rvars(patient_measures, t_offset, stan_data) |> 
    filter(patient_id == 1, tumor_index < 6) |> 
    unnest_rvars()
}

get_tumor_prior_fake_data <- function(res, patient_measures, t_offset, stan_data, ndraws = 24) {
  res |> 
    spread_rvars(rep_tumor_size[tumor_index], ndraws = ndraws) |> # We'll simulate parameters and data for a number simulations. 
    ungroup() |> 
    prepare_tumor_size_rvars(patient_measures, t_offset, stan_data) |> 
    unnest_rvars() |>
    # group_split(.draw, .iteration, .chain)
    nest(fake_data = !c(.draw, .iteration, .chain))
}

sim_tumor_fake_data <- function(data_row, stan_data, tumor_model, output_dir) {
  fit <- stan_data |> 
    list_modify(
      fit_data = TRUE,
      gen_tumor_sizes = FALSE,
      tumor_size = data_row$fake_data[[1]]$tumor_size
    ) |> 
    tumor_model$sample(
      iter_warmup = 200, iter_sampling = 200,
      parallel_chains = 4,
      output_basename = str_c("fake_tumor_", data_row$.draw),
      output_dir = output_dir,
      init = \(chain_id) lst(
        pop_tumor_gp_rho = invgamma::rinvgamma(1, 7.3, 7.5),
        tumor_mean = rnorm(1, 2.8, 0.1),
        tumor_sd = abs(rnorm(1, 0, 1.25))
      ),
      max_treedepth = 15
    )
  
  # fit$save_outputfiles(getOption("cmdstanr_output_dir"), "fake_tumor")
  # return(fit)
  
  data_row |> add_column(fake_fit = list(fit))
}

get_tumor_prior_fake_data_rvar <- function(res_data, tumor_prior_res) {
  res_data |> 
    group_by(.draw) |> 
    reframe(map_dfr(fake_fit, \(f) gather_rvars(f, pop_tumor_gp_intercept, pop_tumor_gp_alpha, pop_tumor_gp_rho))) |> 
    left_join(
      tumor_prior_res |> 
        gather_rvars(pop_tumor_gp_intercept, pop_tumor_gp_alpha, pop_tumor_gp_rho) |> 
        unnest_rvars(), by = c(".draw", ".variable"), suffix = c("", "_true")) |> 
    mutate(.value = .value - .value_true)  
}

get_tumor_analysis_stan_data <- function(treated_tumor_analysis_data, tumor_priors) {
  tumor_stan_data <- lst(
    fit_tumor_data = TRUE,
    gen_tumor_sizes = TRUE,
    predict_missing_sizes = TRUE, 
    multilevel_patient = TRUE,
    multilevel_tumor = FALSE,
    
    !!!prepare_tumor_stan_data(treated_tumor_analysis_data),
    !!!tumor_priors
  )
}

get_all_tumor_size <- function(res, treated_tumor_analysis_data, n_full_measures, .width = 0.8) {
  res |> 
    spread_rvars(all_tumor_size[idx]) |> 
    mutate(
      usubjid = rep(with(treated_tumor_analysis_data, rep(usubjid, n_tumors)), n_full_measures),
      trlnkid = rep(treated_tumor_analysis_data |> unnest(patient_tumors) |> pull(trlnkid), n_full_measures),
    ) |> 
    left_join(
      treated_tumor_analysis_data |> 
        unnest(patient_tumors) |> 
        select(usubjid, trlnkid, trial, matches("(patient_)?(min|max)_t")), 
      by = c("usubjid", "trlnkid")
    ) |> 
    group_by(usubjid, trlnkid) |> 
    mutate(week = first(min_t):first(patient_max_t)) |> 
    ungroup() |> 
    median_qi(all_tumor_size, .width = .width)
}

get_early_fake_tumors <- function(fake_tumor_data) {
  fake_tumor_data |> 
    group_by(patient_id) |> 
    filter(t %in% c(max(keep(t, \(x) x <= 0)), min(keep(t, \(x) x > 0)))) 
}

get_pfs_prior_stan_data <- function(fake_tumor_data, early_fake_tumors, stan_data, n_patients, max_pfs, pfs_priors) {
  stan_data |> 
    list_modify(
      fit_data = FALSE,
      gen_pfs = TRUE,
      gen_interval_censored = FALSE,
      tumor_hazard_type = 1,
      ignore_interval_censoring = FALSE,
      add_trial_level = TRUE,
      add_tumor_location_level = FALSE,
      fit_post_2nd_meaure_only = FALSE, 
      
      pfs = rep(max_pfs, n_patients),
      death_week = rep(0, n_patients),
      right_censored = rep(FALSE, n_patients),
      tumor_size = fake_tumor_data$tumor_size,
      
      grid_tumors = array(NA, dim = 0),
      n_grid_tumors = 0,
      
      tumor_grid_range = with(early_fake_tumors, seq(min(tumor_size), max(tumor_size), 1)), 
      
      !!!pfs_priors,
      
      # !!!tumor_stan_data
      use_tumor_model = FALSE
    ) %>% 
    list_assign(n_tumor_grid_range = length(.$tumor_grid_range))
}

get_max_t <- function(pfs_prior_stan_data) {
  with(pfs_prior_stan_data, list_measures(t_measure, n_measures, n_patient_tumors)) |> 
    map(unlist) |> 
    map_int(max)
}

get_pfs_params <- function(res) {
  res |> 
    spread_rvars(tumor_stim_trial_intercept[trial],
                 tumor_stim_trial_coef[trial, t],
                 log_lambda_gp_alpha,
                 log_lambda_gp_rho,
                 log_lambda_gp_intercept
              ) |>
  pivot_wider(names_from = t, values_from = tumor_stim_trial_coef, names_prefix = "tumor_stim_coef_")
}

get_lambda_rvar <- function(res) {
  res |> 
    spread_rvars(log_lambda[t], base_pf_cond_prob[t], base_survival[t], one_tumor_survival[t], log_trial_lambda[trial, t]) |>
    mutate(lambda = exp(log_lambda), trial_lambda = exp(log_trial_lambda)) |>
    nest(trial_data = contains("trial"))
}

get_expected_pfs_rvar <- function(res) {
  res |> 
    gather_rvars(base_cond_expected_pfs, one_tumor_cond_expected_pfs, base_cond_median_pfs, one_tumor_cond_median_pfs) |>
    mutate(
      stat_type = str_extract(.variable, "expected|median"),
      tumor_size = str_extract(.variable, "base|one_tumor")
    )
}

get_disease_progress_prob <- function(res, stan_data, max_t, n = 2) {
  res |> 
    spread_rvars(disease_progress_prob[i]) |>
    mutate(
      patient_index = rep(seq(stan_data$n_patients), each = max(max_t)), 
      trial = rep(stan_data$patient_trial, each = max(max_t)), 
      t = rep(1:max(max_t), stan_data$n_patients),
    ) |>
    nest(.by = c(trial, patient_index)) |> 
    group_by(trial) |> 
    sample_n(2) |> 
    ungroup() |> 
    unnest(data) |> 
    unnest_rvars()
}

fit_by_trial <- function(trial, stan_data, pfs_model) {
  pfs_model$sample(stan_data, iter_warmup = 300, iter_sampling = 300, 
                  init = create_pfs_initializer(stan_data), parallel_chains = 4,
                  output_dir = file.path(tmp_dir, "fit"), output_basename = str_glue("two_tumor_trial_pfs_{trial}"))
}

bind_pooling_res <- function() {
  map2(two_tumor_trial_loo, split_two_tumor_loo, 
       \(sep_loo, partial_pooled_loo) loo_compare(lst("Fully Separated" = sep_loo, "Partially Pooled" = partial_pooled_loo))) |> 
    map_dfr(\(l) as_tibble(l, rownames = "model"), .id = "trial") |> 
    filter(fct_match(model, "Fully Separated")) |> 
    ggplot(aes(y = trial)) +
    geom_col(aes(x = elpd_diff), fill = AZ_platinum, alpha =0.5, position = "dodge") +
    geom_errorbar(aes(xmin = elpd_diff - se_diff, xmax = elpd_diff + se_diff), width = 0.25) +
    labs(x = "Difference in log-probability score", y = "Trial") +
    NULL
}

get_tumor_design_matrix <- function(stan_data, pfs_functions) { 
  with(stan_data, {
    n_screening_t <- pfs_functions$calc_n_screening_t(n_patient_tumors, n_measures, t_measure)
      
    pfs_functions$standardize_nonzero_tumor_sizes(tumor_size) |>
      pluck(3) |> 
      pfs_functions$prepare_early_tumors_design_matrix(n_patient_tumors, n_measures, t_measure, n_screening_t, 2) |> 
      magrittr::extract2(1)
  })
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

