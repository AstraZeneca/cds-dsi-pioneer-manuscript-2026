gen_patient_interval_properties <- function(log_lambda, tumor_intercept, tumor_coef, settings) {
  standardized <- pfs_model$functions$standardize_nonzero_tumor_sizes(settings$tumor_size)[[3]]
  
  t_measure_list <- with(settings, list_measures(t_measure, n_measures, n_patient_tumors))
  n_screening_t <- with(settings, pfs_model$functions$calc_n_screening_t(n_patient_tumors, n_measures, t_measure)) 
  
  with(settings, pfs_model$functions$prepare_early_tumors_design_matrix(standardized, n_patient_tumors, n_measures, n_screening_t)) |>  
    as_tibble() |> 
    set_names(c("tumor_size_1", "tumor_size_2")) |> 
    mutate(patient_id = rep(1:settings$n_patients, settings$n_patient_tumors)) |> 
    group_by(patient_id) |>
    summarize(tumor_covar = list(cbind(tumor_size_1, tumor_size_2))) |>
    rowwise() |> 
    reframe(
      patient_id, 
      progress_prob = pfs_model$functions$calculate_progress_linear_prob(settings$n_patient_tumors, log_lambda, tumor_intercept, tumor_coef, tumor_covar),
      hazard = pfs_model$functions$calculate_linear_hazard(settings$n_patient_tumors, log_lambda, tumor_intercept, tumor_coef, tumor_covar),
      survival = cumprod(progress_prob)
    )  
}

read_entimice_data <- function(idap, dataset, data_type = c("sdtm", "adam"), team_dir = "/wscratch/ewfteams/dpo0083") {
  read_rds(file.path(team_dir, idap, arg_match(data_type), "prod", "data", str_c(dataset, ".rds"))) |> 
    rename_with(str_to_lower) |> 
    mutate(across(ends_with("fl"), \(fl) fct_expand(fl, c("Y", "N")) |>  fct_match("Y")))
}

km_to_tibble <- function(trt_data, key, pfs_var) { 
  with(
    prepare_pfs_stan_data(trt_data, tumor_priors, pfs_priors, pfs_var = pfs_var), {
      interval_censored <- pfs_model$functions$identify_censoring(pfs, death_week, n_patient_tumors, n_measures, t_measure)[[1]]
      
      map_dfr(list(lb = pfs, ub = pfs + interval_censored), function(s) {
        pfs_model$functions$estimate_kaplan_meier(s, right_censored, max(s)) |>
          set_names(c("s", "n", "c", "e")) |>
          as_tibble() |> 
          mutate(t = seq(0, n() - 1))
      }, .id = "btype")
    }) |> 
    # pivot_wider(values_from = s:e, names_from = btype) |> 
    bind_cols(key)
}

list_measures <- function(measures, n_measures, n_patient_tumors) { 
  split(measures, rep(seq_along(n_measures), n_measures)) |>  
    split(rep(seq_along(n_patient_tumors), n_patient_tumors))
}

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

gen_fake_pfs_data <- function(patient_interval_data, settings) { 
  fake_data <- patient_interval_data |>
    nest(prob = !patient_id) |>  
    mutate(t_measure = with(settings, list_measures(t_measure, n_measures, n_patient_tumors))) %>% 
    mutate(
      map2(.$prob, .$t_measure, \(pd, t) pfs_model$functions$pfs_rng(pd$progress_prob, discard(first(t), \(x) x <= 0))) |> 
        list_transpose() |> 
        set_names(c("interval_censored", "right_censored", "pfs", "actual_pfs")) |> 
        `!!!`(),
    ) %>%
    mutate(
      pfs_model$functions$identify_censoring(.$pfs, rep_along(.$pfs, FALSE), settings$n_patient_tumors, settings$n_measures, settings$t_measure) |> 
        set_names(c("stan_interval_censored", "stan_right_censored")) |> 
        `!!!`()
    )
  
  assertthat::assert_that(with(fake_data, all(stan_interval_censored == interval_censored)))
  assertthat::assert_that(with(fake_data, all(stan_right_censored == right_censored)))
  
  return(fake_data)
}

create_pfs_initializer <- function(stan_data) {
  function(chain_id) { 
    n_covar <- with(stan_data, if_else(tumor_hazard_type > 0 && tumor_hazard_type != 4, 
                                       if_else(tumor_hazard_type < 3, 
                                               tumor_hazard_type + 1, 
                                               if_else(tumor_hazard_type == 5, 5, 1)), 
                                       0))
    
    init_vals <- lst(
      tumor_stim_intercept = abs(rnorm(1, 0, stan_data$tumor_stim_intercept_sd)),
      tumor_stim_coef = abs(rnorm(n_covar, 0, stan_data$tumor_stim_coef_sd[1:n_covar])), 
    )
    
    if (stan_data$add_trial_level) {
      init_vals <- init_vals |> 
        list_assign(
          raw_log_lambda_gp_trial_intercept = rnorm(stan_data$n_trials),
          raw_tumor_stim_trial_coef_mult = map(seq(stan_data$n_trials), \(...) rnorm(n_covar + 1)),
          log_lambda_gp_trial_intercept_sd = abs(rnorm(1, 0, stan_data$log_lambda_gp_trial_intercept_sd_sd)),
          tumor_stim_trial_coef_mult_sd = abs(rnorm(n_covar + 1, 0, stan_data$tumor_stim_trial_coef_sd_sd)),
        )
    }
    
    return(init_vals)
  }
}

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
      init = if (settings$tumor_hazard_type > 0 && !no_init) create_pfs_initializer(settings),
      ...
    )
}

fit_simulations <- function(n, patient_interval_data, settings, ignore_interval_censoring = FALSE, fit_basename = NULL, tmp_dir = here("temp")) {
  fake_data_sim_with_ic <- tibble(sim_id = seq(n)) |> 
    rowwise() |> 
    mutate(sim_data = list(gen_fake_pfs_data(patient_interval_data, settings))) |> 
    ungroup() |> 
    transmute(
      sim_id,
      sim_data,
      sim_fit = furrr::future_map2(.progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
        sim_id, sim_data, 
        \(sid, sdata) fit_sim_data(
          settings, sdata, ignore_interval_censoring = ignore_interval_censoring, 
          output_basename = if (!is_null(fit_basename)) str_c(fit_basename, sid, sep = "_"), 
          output_dir = if (!is_null(fit_basename)) file.path(tmp_dir, "fit"))
      ), 
    ) 
}

prepare_tumor_stan_data <- function(analysis_data) {
  lst(
    n_patients = nrow(analysis_data),
    n_trials = n_distinct(analysis_data$trial),
    patient_trial = analysis_data$trial,
    n_patient_tumors = analysis_data$n_tumors,
    n_measures = analysis_data$n_measures |> unlist(),
    t_measure = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$week) |> unlist(),
    tumor_size = unnest(analysis_data, patient_tumors) |> pull(tumor_history) |> map(\(h) h$mmdiam / 10) |> unlist(),
  )
}

prepare_tumor_size_rvars <- function(rv, settings) {
  rv |> 
    ungroup() |> 
    mutate( # A bunch of acrobatics to get the IDs, indices, and intervals right.
      patient_id = rep(1:n_patients, with(settings, n_patient_tumors * patient_measures)),
      tumor_index = ((tumor_index - 1) %/% patient_measures) + 1,
      t = rep(0:(patient_measures - 1), sum(settings$n_patient_tumors)) - t_offset
    ) |>
    mutate(tumor_index = tumor_index - min(tumor_index) + 1, .by = patient_id) |>
    rename(tumor_size = rep_tumor_size)
}

prepare_pfs_stan_data <- function(analysis_data, .tumor_priors, .pfs_priors, ..., pfs_var = pfs) {
    tumor_stan_data <- prepare_tumor_stan_data(analysis_data)
    pfs_data <- select(analysis_data, pfs = {{ pfs_var }}, death_week, right_censored, interval_censored) |> 
      mutate(death_week = if_else(right_censored, 0, death_week)) # Death week is irrelevant if the data is censored
  
    lst(
    fit_data = TRUE,
    use_tumor_model = FALSE,
    gen_pfs = TRUE,
    gen_interval_censored = FALSE,
    ignore_interval_censoring = FALSE,
    add_trial_level = FALSE,
    
    fit_tumor_data = FALSE,
    gen_tumor_sizes = FALSE,
    predict_missing_sizes = FALSE, 
    multilevel_patient = FALSE,
    multilevel_tumor = FALSE,
    
    !!!tumor_stan_data,
    !!!pfs_data,
    
    !!!.pfs_priors,
    !!!.tumor_priors,
    
    ...
  ) |> 
    list_modify(fit_tumor_data = FALSE)
}

run_sbc_sims <- function(
  pfs_model, stan_data, num_sim, output_name,  output_dir = file.path(tmp_dir, "fit"), 
  ignore_interval_censoring = FALSE,keep_fit = FALSE, gen_pfs = FALSE,
  reuse_data = NULL, drop_measures = NULL, keep_only = FALSE, ... 
) {
  spread_param_rvars <- function(f, ...) { 
    f |> 
      spread_rvars(
        tumor_stim_intercept, tumor_stim_coef[t],
        log_lambda_gp_intercept, log_lambda_gp_alpha, log_lambda_gp_rho,
        base_cond_expected_pfs, one_tumor_cond_expected_pfs, base_cond_median_pfs, one_tumor_cond_median_pfs,
        ...
      ) |> 
      ungroup() |> 
      pivot_wider(names_from = t, values_from = tumor_stim_coef, names_prefix = "tumor_stim_coef_")
  }
  
  stan_data <- stan_data |> 
    drop_missing_measures(drop_measures, keep_only)
  
  lstm <- with(stan_data, list_measures(t_measure, n_measures, n_patient_tumors))
 
  all_sim_data <- if (is_null(reuse_data)) { 
    pfs_res <- pfs_model$sample(data = stan_data, refresh = 0, parallel_chains = 4)
    
    pfs_res |> 
      spread_rvars(rep_pfs[patient_index], rep_right_censored[patient_index]) |> 
      unnest_rvars() |> 
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
  } else {
    reuse_data |> select(.draw, sim_data, true)
  }
  
  all_sim_data |> 
    mutate(# For each simulation dataset fit the model and extract the posteriors for each of the model parameters 
      sim_data = map(
        sim_data, 
        \(d) mutate(d, 
                    map2_dfr(pfs, lstm, function(s, m) { 
                      m_union <- reduce(m, \(a, n) union(a, n)) 
                      list(pfs = m_union |> discard(\(t) t > s) |> max(), right_censored = s >= max(m_union))
                    })
                    # pfs = map2_int(pfs, lstm, \(s, m) reduce(m, \(a, n) union(a, n)) |> discard(\(t) t > s) |> max()),
                    # right_censored = pmap_lgl(lst(c = right_censored, p = pfs, m = lstm), \(c, p, m) c || )
        )),
                              
      furrr::future_map2_dfr(.draw, sim_data, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE),
      # map2_dfr(.draw, sim_data,
        function(sim_id, d, output_dir, output_name, ignore_interval_censoring) { 
          fit <- fit_sim_data(
            stan_data,
            d, 
            gen_pfs = gen_pfs, 
            thin = 4, # We need thinning when doing SBC using MCMC to break the correlation between samples.
            output_basename = str_glue("{output_name}_{sim_id}"),
            output_dir = output_dir, 
            ignore_interval_censoring = ignore_interval_censoring,
            ...
          ) 
          
          res <- fit |> # Get posterior draws from simulation fit. 
            spread_param_rvars(ndraws = 1000) |> 
            pack(est = everything())
          
          if (keep_fit) {
            res <- res |> 
              mutate(sim_fit = list(fit))
          }
          
          return(res)
        },
        output_dir = output_dir, output_name = output_name, ignore_interval_censoring = ignore_interval_censoring
      ),
    ) 
}

plot_km <- function(data_list, analysis_data, facet_arm = FALSE) {
  plot_obj <- data_list |> 
    map(\(r) recover_types(r, select(analysis_data, trial))) |> 
    map_dfr(\(r) spread_rvars(r, trial_km_est[trial, t]), .id = "arm") |> 
    ggplot() +
    labs(title = "Kaplan-Meier estimate", subtitle = "Treated arm", x = "t", y = latex2exp::TeX("$S(t)$")) +
    NULL
  
  if (!facet_arm) {
    plot_obj <- plot_obj + facet_wrap(vars(trial)) 
  } else {
    plot_obj <- plot_obj + facet_grid(vars(arm), vars(trial)) 
  }
  
  if (length(data_list) > 1 && !facet_arm) {
    plot_obj +
      stat_lineribbon(aes(x = t - 1, ydist = trial_km_est, fill = arm, color = arm), step = TRUE, alpha = 0.25, .width = 0.8) 
  } else {
    plot_obj +
      stat_lineribbon(aes(x = t - 1, ydist = trial_km_est), fill = "black", step = TRUE, alpha = 0.25, .width = 0.8) 
  }
}

plot_pfs_hist_posterior <- function(data_list, stan_data, hist_breaks = seq(10, 150, 10)) {
  # This function is used to generate a histogram of time-to-events for a single draw
  sample_hist <- function(pred, breaks) {
    # hist() is a base R function to generate histograms from data and provided breaks.
    hist(pmin(pred, max(breaks)), breaks = c(0, breaks), plot = FALSE)$count
  }
  
  # This function is used to allow us to generate a distribution of histograms
  rvar_sample_hist <- posterior::rfun(sample_hist)
  
  data_list |> 
    map_dfr(\(f) spread_rvars(f, rep_pfs[i], rep_right_censored[i]) |> mutate(trial = stan_data$patient_trial), .id = "arm") |> 
    group_by(trial) |> 
    reframe(t = hist_breaks, bindist = rvar_sample_hist(rep_pfs, hist_breaks)) |> 
    filter(t < max(t)) %>% 
    bind_rows(
      group_by(., trial) %>%
        filter(t %in% range(t)) |> 
        mutate(t = c(0, max(t) + min(hist_breaks)))
    ) |> 
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = bindist, color = "Posterior"),
                    fill = "black", alpha = 0.125, linewidth = 2,
                    step = "mid",
                    .width = c(0.5, 0.8), show.legend = FALSE) + 
    scale_x_continuous("t", breaks = seq(0, max(hist_breaks) + min(hist_breaks), 20)) +
    scale_color_discrete("") +
    facet_wrap(vars(trial)) +
    NULL
}

plot_base_hazard <- function(data_list, analysis_data) {
  data_list |> 
    map(\(f) recover_types(f, select(analysis_data, trial))) |> 
    map_dfr(\(f) spread_rvars(f, log_trial_lambda[trial, t]), .id = "fit_type") |> 
    mutate(trial_lambda = exp(log_trial_lambda)) |> 
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = trial_lambda, color = fit_type, fill = fit_type, alpha = fit_type), step = TRUE, .width = 0.8) +
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

plot_loss <- function(data, binwidth) {
  data |> 
    ggplot(aes(value)) +
    stat_histinterval(
      aes(fill = model, color = model), point_interval = mean_qi, .width = c(0.5, 0.8), alpha = 0.5, breaks = ggdist::breaks_fixed(width = binwidth)
    ) +
    # scale_color_viridis_d(option = "E", aesthetics = c("color", "fill")) +
    labs(x = "Loss", y = "") +
    facet_grid(vars(trial), vars(loss_type), scales = "free", margins = "trial") +
    theme(legend.position = "top", strip.text.y = element_text(angle = 0)) +
    NULL
}
