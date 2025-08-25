# nolint start: object_usage_linter

prepare_pdl1_and_trial_info <- function(res_data) {
  res_data |>
    filter(!fct_match(variable, "first_liners") | !fct_match(cond_group_name, "no")) |> 
    mutate(
      trial = coalesce(trial, "sclc"),
      across(c(cond_group_name, variable), \(l) coalesce(l, "all")),
      variable = fct_collapse(variable, "all" = c("all", "pdl1"), "pdl1_naive" = c("pdl1_naive", "first_liners")),
      cond_group_name = fct_collapse(cond_group_name, "All" = c("all", "yes"), "PDL 1 Low" = "low", "PDL1 High" = "hi")
    )
}

#' Plot outcome (e.g., ORR or Median PFS) by PDL1 status and trial for SCLC trial
#'
#' @param data Data frame with columns: outcome, cond_group_name, variable, trial, fit_type, etc.
#' @param outcome Unquoted column name for the outcome to plot (e.g., orr, median_pfs)
#' @param ... Additional arguments passed to ggplot2::stat_pointinterval
#' @return A ggplot object
plot_outcome_by_pdl1_and_trial <- function(res_data, outcome, ...) {
  res_data |> 
    prepare_pdl1_and_trial_info() |>
    filter(fct_match(variable, c("all", "pdl1_naive")), fct_match(trial, "sclc")) |>  
    ggplot() +
    stat_pointinterval(aes(xdist = {{ outcome }}, y = cond_group_name, color = fit_type), position = "dodge", .width = c(0.5, 0.9), ...) +
    scale_color_discrete("", type = AZ_palette, label = str_to_title) +
    facet_grid(vars(variable), vars(trial), scales = "free", space = "free", 
               labeller = labeller(trial = str_to_upper, variable = c("all" = "All", "pdl1_naive" = "First Line"))) +
    NULL
}

plot_baseline_hazard <- function(res_data, lambda_var, ...) {
  ggplot(res_data, aes(t)) +
    geom_line(
      aes(y = {{ lambda_var }}, color = fit_type, group = str_c(fit_type, .draw)),
      alpha = 0.25,
      data = \(d) unnest_rvars(d) |>
        ungroup() %>%
        semi_join(
          distinct(., fit_type, ..., .draw) |>
            group_by(fit_type, ...) |>
            slice_sample(n = 25),
          by = join_by(fit_type, ..., .draw)
        )
    ) +
    stat_lineribbon(aes(ydist = {{ lambda_var }}, fill = fit_type, color = fit_type), alpha = 0.25, linewidth = 0, .width = 0.8) +
    labs(y = "Baseline Hazard") +
    theme(legend.position = "bottom")
}

plot_crcr_baseline_hazard <- function(res_data, analysis_data = NULL) {
  po <- plot_baseline_hazard(res_data, crcr_trial_lambda, k) +
    facet_grid(vars(trial), vars(k), scales = "free_y")
  
  if (!is_null(analysis_data)) {
    po <- po +
      geom_rug(
        aes(x = confirmed_response_week),
        alpha = 0.5,
        data = analysis_data |>
          filter(!confirmed_response_censored) |>
          mutate(k = if_else(confirmed_response, "Response", "Non-response"))
      )
  }
  
  return(po)
}

plot_pfs_baseline_hazard <- function(res_data, analysis_data) {
  plot_baseline_hazard(res_data, trial_lambda) +
    geom_rug(aes(x = progress_week), alpha = 0.5, data = analysis_data |> filter(!right_censored)) +
    labs(y = "Baseline Hazard") +
    facet_wrap(vars(trial)) +
    theme(legend.position = "bottom") 
}

plot_unclassified_survival <- function(res_data, analysis_data, conf_resp_hb) {
  ggplot(res_data, aes(t)) +
    stat_lineribbon(
      aes(ydist = bindist, color = fit_type, fill = fit_type, alpha = fit_type),
      linewidth = 1, step = "hv", .width = c(0.5, 0.8)
    ) +
    geom_step(
      aes(y = count, linetype = "Observed"), direction = "vh", linewidth = 0.5,
      data = \(d) analysis_data |>
        group_by(trial) |> 
        reframe(t = conf_resp_hb[-length(conf_resp_hb)], count = sample_hist(confirmed_response_week, conf_resp_hb))
    ) +
    scale_alpha_manual("", values = c(prior = 0.125, posterior = 0.25)) +
    scale_linetype_manual("", values = c(Observed = "dotted")) +
    labs(y = "Count") +
    theme(legend.position = "bottom") +
    guides(alpha = "none") +
    NULL 
}

plot_cif <- function(res_data, obs_cif_data) {
  ggplot(res_data) +
    stat_lineribbon(aes(x = t, ydist = trial_cif, fill = fit_type), linewidth = 0, alpha = 0.25, .width = c(0.5, 0.8)) +
    geom_step(aes(x = time, y = estimate, linetype = "Observed"), direction = "vh", data = \(d) semi_join(obs_cif_data, d, by = "trial")) +
    scale_linetype_manual("", values = c(Observed = "dashed")) +
    labs(y = "CIF") +
    NULL
}

plot_hazard_ratio <- function(res_data) {
  ggplot(res_data, aes(t)) +
    stat_lineribbon(aes(ydist = bindist, color = fit_type, fill = fit_type, alpha = fit_type),
                    linewidth = 1,
                    step = "hv",
                    .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    scale_alpha_manual("", values = c(prior = 0.125, posterior = 0.25)) +
    labs(x = "Hazard Ratio", y = "Density") +
    guides(alpha = "none") +
    theme(legend.position = "bottom", strip.text.y.left = element_text(angle = 0), strip.placement = "outside") +
    NULL 
}

plot_crcr_covar_coef <- function(res_data) {
  res_data |> 
    filter(fct_match(.variable, "crcr_covar_trial_coef")) |> 
    select(!.value) |> 
    pivot_wider(names_from = k, values_from = .exp_value) |> 
    mutate(tumor = str_detect(covar, "tumor sizes"), hazard_ratio = Response / `Non-response`) |> 
    ggplot(aes(y = covar)) +
    stat_pointinterval(aes(xdist = hazard_ratio, color = fit_type), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    labs(x = "Ratio of Competing Risk Hazards",  y = "Parameter", caption = "Showing the posterior median, 50% CI, and 80% CI.") +
    coord_cartesian(clip = "off") +
    theme(legend.position = "bottom", strip.text.y = element_blank()) +
    NULL 
}

plot_pfs_covar_coef <- function(res_data) {
  res_data |> 
    filter(fct_match(.variable, "covar_trial_coef")) |> 
    ggplot(aes(y = covar)) +
    stat_pointinterval(aes(xdist = .exp_value, color = fit_type), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    labs(x = "Exponential of Parameter",  y = "Parameter", caption = "Showing the posterior median, 50% CI, and 80% CI.") +
    theme(legend.position = "bottom", strip.text.y = element_blank()) +
    NULL
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

plot_km <- function(res_data, obs_km_data, km_est, analysis_data = NULL, group = fit_type, alpha_group = fit_type, color_group = fit_type, linewidth = 0, ...) {
  pobj <- ggplot(res_data, aes(x = t - 1)) +
    stat_lineribbon(
      aes(
      dist = {{ km_est }},
      fill = {{ group }},
      alpha = {{ alpha_group }},
      group = {{ group }},
      ),
      linewidth = linewidth, .width = 0.8, ...
    ) +
    # stat_lineribbon(
    #   aes(
      # dist = {{ km_est }},
    #   color = {{ color_group }}
    #   ),
    #   .width = 0, # only median line
    #   linewidth = linewidth,
    #   alpha = 1,
    #   fill = NA,
    #   show.legend = FALSE
    # ) +
    labs(y = "Survival Probability") +
    guides(alpha = "none") + 
    theme(legend.position = "bottom")
  
  if (!is_null(obs_km_data)) {
    pobj <- pobj + 
      geom_step(aes(x = t, y = s, group = btype, color = btype), linewidth = 0.5, alpha = 0.5, data = \(d) semi_join(obs_km_data, d, by = "trial"))
    
    if (!is_null(analysis_data)) {
      pobj <- pobj +
        geom_point(aes(x = t, y = s, color = btype, shape = "censored"), size = 2, alpha = 0.7, 
                   data = \(d) semi_join(obs_km_data, d, by = "trial") |> 
                     inner_join(analysis_data |> filter(right_censored) |> select(pfs), by = c("t" = "pfs"), relationship = "many-to-many")) +
        geom_point(aes(x = t, y = s, color = btype, shape = "death"), size = 2, alpha = 0.7, 
                   data = \(d) semi_join(obs_km_data, d, by = "trial") |> 
                     inner_join(analysis_data |> filter(!right_censored, !progression_before_death) |> select(pfs), by = c("t" = "pfs"), relationship = "many-to-many")) +
        scale_shape_manual("", values = c(censored = "|", death = "o"), labels = c(censored = "Right Censored", death = "Death before PD"))
    }
  }
  
  return(pobj)
}

plot_gng <- function(res_data, outcome, lrv_tv, model_type_names) {
  res_data |> 
    ggplot(aes(y = data_cut)) +
    stat_interval(aes(xdist = {{ outcome }}, color = impute_type, color_ramp = after_stat(level)),  position = "dodge", .width = c(0.6, 0.8)) +
    geom_vline(xintercept = lrv_tv, linetype = "dashed") +
    scale_x_continuous("", sec.axis = sec_axis(identity, breaks = lrv_tv, labels = c("LRV", "TV"))) +
    scale_y_discrete("", labels = str_to_title) +
    scale_color_discrete("Sample", type = AZ_palette, labels = \(l) str_replace(l, "_", " ") |> str_to_title()) +
    scale_color_ramp_discrete(name = "Credible Intervals") +
    facet_grid(vars(model_type), switch = "y", labeller = labeller(model_type = model_type_names)) +
    theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
}

plot_simple_gng <- function(res_data, outcome, color_col, lrv_tv, model_type_names, y_col = model_type, outcome_desc = "", decision_prob = c(0.2, 0.9)) {
  res_data |> 
    ggplot(aes(y = {{ y_col }})) +
    stat_interval(
      aes(xdist = {{ outcome }}, color = {{ color_col }}, color_ramp = after_stat(level)), 
      position = "dodge", .width = c(0.6, 0.8)
    ) +
    geom_segment(
      aes(x = decision_interval[, 1], xend = decision_interval[, 2], group = {{ color_col }}), 
      arrow = grid::arrow(ends = "both", length = unit(1.5, "mm")), 
      position = position_dodge(width = 1),
      data = \(d) mutate(d, decision_interval = t(quantile({{ outcome }}, decision_prob)))
    ) +
    geom_vline(xintercept = lrv_tv, linetype = "dashed") +
    scale_x_continuous(outcome_desc, sec.axis = sec_axis(identity, breaks = lrv_tv, labels = c("LRV", "TV"))) +
    scale_y_discrete("", labels = model_type_names) +
    scale_color_discrete("Sample", type = AZ_palette, labels = \(l) str_replace(l, "_", " ") |> str_to_title()) +
    scale_color_ramp_discrete(name = "Credible Intervals", range = c(0.25, 0.5)) +
    NULL
}

add_dco_to_plot <- function(plot, label_offset_x = - days(40), label_offset_y = 2) {
  plot +
    geom_vline(xintercept = c(lubridate::ymd("2024-04-22"), lubridate::ymd("2024-07-31"), lubridate::ymd("2024-10-31")), linetype = "dashed", color = AZ_platinum) +
    annotate("text", x = lubridate::ymd("2024-04-22") + label_offset_x, y = label_offset_y, label = "DCO 0") +
    annotate("text", x = lubridate::ymd("2024-07-31") + label_offset_x, y = label_offset_y, label = "DCO 1") +
    annotate("text", x = lubridate::ymd("2024-10-31") + label_offset_x, y = label_offset_y, label = "DCO 2") 
}

plot_patient_timelines <- function(analysis_data) {
  plot <- analysis_data |> 
    unnest(patient_tumors) |>
    unnest(tumor_history, names_sep = "_") |> 
    distinct(usubjid, trtsdt, right_censored, day = tumor_history_day, visit_date = tumor_history_adt) |> 
    nest(visits = c(day, visit_date)) |> 
    mutate(map_dfr(visits, \(v) summarize(v, first_visit = min(visit_date), last_visit = max(visit_date)))) |> 
    mutate(usubjid = fct_reorder(usubjid, first_visit)) |> 
    ggplot(aes(y = usubjid)) +
    geom_segment(aes(x = first_visit, xend = last_visit, yend = usubjid)) +
    geom_point(aes(x = visit_date, shape = "visit"), size = 2, data = \(d) unnest(d, visits) |> filter(visit_date < last_visit)) +
    geom_point(aes(x = trtsdt, shape = "treat"), size = 2) +
    geom_point(aes(x = last_visit, color = right_censored, shape = "last"), size = 2) +
    labs(x = "Calendar Time", y = "Patients") +
    scale_color_discrete("", label = c("FALSE" = "Progression", "TRUE" = "Censored"), type = AZ_palette) +
    scale_shape_manual(
      "", values = c("visit" = 124, "treat" = 5, "last" = 19), labels = c("visit" = "Visit", "treat" = "Treatment Start", "last" = "Last Visit")
    ) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "inside", legend.position.inside = c(0.25, 0.8))
  
  add_dco_to_plot(plot) 
}

get_pfs_conf_resp_marginal_exit_prob <- function(res) {
  res |>
    ungroup() |> 
    transmute(trial, prob_rvars = map(fit, \(f) spread_rvars(f, marginal_exit_prob[i, k, t]))) |> 
    unnest(prob_rvars)
}

get_trial_sim_crcr_pfs <- function(res, stan_data) {
  spread_rvars(res, sim_pfs[i], sim_censored[i], forecast_pfs[i], forecast_censored[i]) |>
    mutate(usubjid = stan_data$patient)
}

get_sim_pfs_conf_resp <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      sim_pfs = map2(
        fit, analysis_data, 
        \(f, d) spread_rvars(f, sim_pfs[i], sim_censored[i]) |>
          left_join(transmute(d, i = seq(n()), pfs, right_censored), by = "i", relationship = "one-to-one")
      )
    ) |> 
    unnest(sim_pfs)
}

get_pfs_conf_resp_km_est <- function(res) {
  res |> 
    select(trial, fit) |> 
    deframe() |> 
    map_dfr(\(r) spread_rvars(r, km_est[t]), .id = "trial") 
}

get_all_pfs_conf_resp_km_est <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
  spread_rvars(res, trial_km_est[trial, t], forecast_trial_km_est[trial, t])
}

get_median_pfs_conf_resp <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, sim_median_pfs))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_median_pfs_conf_resp <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
   spread_rvars(res, sim_trial_median_pfs[trial], forecast_trial_median_pfs[trial]) 
}

get_pfs_n <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
   spread_rvars(res, sim_trial_pfs6[trial], sim_trial_pfs9[trial], forecast_trial_pfs6[trial], forecast_trial_pfs9[trial]) 
}

get_pfs_conf_resp_log_hazard_ratio <- function(res) {
   res |> 
     rowwise() |> 
     transmute(trial, rv = list(spread_rvars(fit, time_invariant_log_hazard_ratio[i, k]) |> 
                                  mutate(time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio)))) |> 
     ungroup() |> 
     unnest(rv)
}

get_all_pfs_conf_resp_hazard_ratio <- function(res, stan_data) {
  spread_rvars(res, time_invariant_log_hazard_ratio[i, k]) |>
    mutate(
      time_invariant_hazard_ratio = exp(time_invariant_log_hazard_ratio),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response")) 
    ) |> 
    left_join(as_tibble(stan_data["patient_trial"]) |> mutate(i = seq(n())), by = "i", relationship = "many-to-one") |> 
    rename(trial = patient_trial)  
}

get_all_pfs_conf_resp_hazard_ratio_bindist <- function(res, stan_data, hb) {
  get_all_pfs_conf_resp_hazard_ratio(res, stan_data) |> 
    group_by(trial, k) |> 
    reframe(t = hb[-length(hb)], bindist = rvar_sample_hist(time_invariant_hazard_ratio, hb))  
}

get_pfs_conf_resp_bootstrap_variables <- function(res, bs_sample1, bs_sample2, rates_name, ..., summarize = TRUE) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(
        spread_rvars(fit, ...) |>  
          mutate(
            n_bs_sample = {{ bs_sample1 }} + {{ bs_sample2 }},
            across(ends_with("predicted"), \(n)  n / n_bs_sample, .names = "{.col}_prop")
          ) %>% { 
            if (summarize) {
              unnest_rvars(.) |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
                point_interval(
                  na.rm = TRUE, # if n_bs_sample is 0, we'll get some NaNs. These are few so we'll bite the bullet and drop them.
                  .width = c(0.5, 0.8)
                ) 
            } else .
          } |> 
          left_join(as_tibble(stan_data[rates_name]) |> mutate(r = seq(n())), by = "r", relationship = "many-to-one")
      )
    ) |> 
    ungroup() |> 
    unnest(rv)  
}

get_pfs_conf_resp_bootstrap_cr_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, "bootstrap_cr_maturity_rates",
    bs_cr_prediction_calendar_week[r], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r],
    bs_cr_median_pfs[r], bs_cr_orr[r],
    n_bs_cr_conf_resp_predicted[r], n_bs_cr_pfs_predicted[r],
    summarize = FALSE
  ) |> 
    mutate(log_bs_cr_median_pfs = log(bs_cr_median_pfs)) |> 
    unnest_rvars() |> # na.rm = TRUE doesn't work in point_interval() if using rvars.  
    point_interval(na.rm = TRUE, .width = c(0.5, 0.8))
}

get_pfs_conf_resp_bootstrap_pfs_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, "bootstrap_pfs_maturity_rates",
    bs_cr_conf_resp_censored_prop[r], bs_pfs_prediction_calendar_week[r], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r],
    bs_pfs_median_pfs[r], bs_pfs_orr[r],
    n_bs_pfs_conf_resp_predicted[r], n_bs_pfs_pfs_predicted[r]
  )
}    

get_fixed_bootstrap_cr_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified, "bootstrap_cr_maturity_rates",
    fixed_bs_cr_median_pfs[r, f], fixed_bs_cr_orr[r, f], n_bs_sample_cr_classified[r], n_bs_sample_cr_unclassified[r],
    n_fixed_bs_cr_conf_resp_predicted[r, f], n_fixed_bs_cr_pfs_predicted[r, f],
    summarize = FALSE
  ) |> 
    select(!c(n_bs_sample, n_bs_sample_cr_classified, n_bs_sample_cr_unclassified))
}

get_fixed_bootstrap_pfs_median_pfs <- function(res) {
  get_pfs_conf_resp_bootstrap_variables(
    res, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving, "bootstrap_pfs_maturity_rates",
    fixed_bs_pfs_median_pfs[r, f], fixed_bs_pfs_orr[r, f], n_bs_sample_pfs_progressed[r], n_bs_sample_pfs_surviving[r],
    n_fixed_bs_pfs_conf_resp_predicted[r, f], n_fixed_bs_pfs_pfs_predicted[r, f],
    summarize = FALSE
  ) |> 
    select(!c(n_bs_sample, n_bs_sample_pfs_progressed, n_bs_sample_pfs_surviving))
}
  
get_sample_maturity_rvar <- function(res) {
  res |>
    ungroup() |> 
    transmute(
      trial,
      rv = map2(
        fit, stan_data, 
        \(f, d) spread_rvars(f, n_sample[l, p], maturity_rate[l, p]) |>
          bind_cols(expand.grid(d[c("lambda", "pred_week")]))
      )
    ) |> 
    unnest(rv)
}

get_all_pfs_crcr_trial_lambda_residual <- function(res) {
  spread_rvars(res, log_trial_lambda_residual[trial, t]) |> 
    mutate(trial_lambda_residual = exp(log_trial_lambda_residual)) |> 
    point_interval(log_trial_lambda_residual, trial_lambda_residual, .width = c(0.5, 0.8))  
}


get_all_pfs_crcr_trial_lambda_residual_draws <- function(res, ndraws = NULL) {
  spread_rvars(res, log_trial_lambda_residual[trial, t]) |> 
    mutate(
      log_trial_lambda_residual = thin_draws(log_trial_lambda_residual),
      trial_lambda_residual = exp(log_trial_lambda_residual)
    ) |> 
    unnest_rvars() |> 
    filter(is_null(ndraws) | (.draw <= ndraws))
}

get_pfs_pred_param <- function(res) {
  res |> 
    rowwise() |> 
    transmute(trial, rv = list(get_all_pfs_pred_param(fit))) |> 
    unnest(rv)  
}

get_cr_median_pfs_draws <- function(res, ndraws = Inf) {
  res |> 
    rowwise() |> 
    transmute(
      trial, 
      rv = list(spread_draws(fit, bs_cr_median_pfs[r]) |>
                  filter(.draw <= ndraws) |> # I use this to make sure all r have the same .draw 
                  left_join(as_tibble(stan_data["bootstrap_cr_maturity_rates"]) |> 
                              mutate(r = seq(n())), 
                            by = "r", relationship = "many-to-one"))
    ) |> 
    unnest(rv)
}

get_all_pfs_crcr_lambda <- function(res, stan_data = NULL) {
  spread_rvars(res, log_trial_lambda[trial, t]) |> 
    mutate(
      trial_lambda = exp(log_trial_lambda), 
      trial = if (!is_null(stan_data)) factor(trial, labels = levels(stan_data$patient_trial))
    )
}

get_all_pfs_crcr_lambda_trial_intercept <- function(res) {
  spread_rvars(res, log_lambda_gp_trial_intercept[trial]) |> 
    mutate(lambda_gp_trial_intercept = exp(log_lambda_gp_trial_intercept))
}

get_crcr_pfs_pred_param <- function(res, stan_data) {
  gather_rvars(res, covar_trial_coef[trial, m], covar_effect[trial, m], tumor_stim_pop_coef[trial, m]) |> 
    mutate(.exp_value = exp(.value)) |> 
    name_coef_indices(m, trial, stan_data)
}

get_powerscaled_variables <- function(res, metadata, stan_data) {
  metadata |> 
    rowwise() |> 
    mutate(
      ps = list(
        if (alpha == 1) res else powerscale(res, alpha = alpha, component = component, variable = c("log_crcr_trial_lambda", "log_trial_lambda"))
      )
    ) |> 
    transmute(
      alpha, component,
      baseline_hazard_rvar = list(
        gather_rvars(ps, log_crcr_trial_lambda[k, trial, t], log_trial_lambda[trial, t]) |> 
          mutate(
            .exp_value = exp(.value),
            trial = factor(trial, labels = levels(stan_data$patient_trial)),
            k = factor(k, levels = 1:2, labels = c("Non-response", "Response")) 
          ) |>  
          point_interval(.value, .exp_value, .width = c(0.5, 0.8))
      ) 
    )
}

plot_coef <- function(d, xvar) {
  ggplot(d) + 
    stat_slab(aes(xdist = {{ xvar }}, color = fit_type), fill = NA, linewidth = 2, show.legend = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 2) +
    scale_color_discrete("", label = str_to_title, type = AZ_palette, aesthetic = c("color", "fill")) +
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 25), legend.text = element_text(size = 20)) +
    coord_cartesian(xlim = c(0, 4)) +
    theme(legend.position = "bottom") +
    NULL
}

get_coef_powerscale_table_data <- function(coef_ps_sense, prior_crcr_coef, crcr_coef, prior_crcr_pfs_coef, crcr_pfs_coef, stan_data, trial_col = NULL) {
  coef_plots <- bind_rows(
    bind_rows(prior = prior_crcr_coef, posterior = crcr_coef, .id = "fit_type"),
    bind_rows(prior = prior_crcr_pfs_coef, posterior = crcr_pfs_coef, .id = "fit_type")
  ) |>
    filter(fct_match(.variable, c("crcr_covar_effect", "covar_effect", "crcr_tumor_stim_pop_coef", "tumor_stim_pop_coef"))) |> 
    nest(coef_data = !c(.variable, m, k)) |> 
    mutate(plot_obj = map(coef_data, \(d) plot_coef(d, .exp_value))) 
  
  coef_ps_sense |> 
    filter(!str_detect(variable, "trial_sd")) |> 
    tidyr::separate_wider_regex(variable, c(var = ".+", r"{\[}", trial = r"{\d+}", ",", m = r"{\d+}", r"{,?}", k = r"{(?:\d+)?}", ".*")) |> 
    mutate(
      across(c(m, k, trial), as.integer),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))
    ) |> 
    left_join(coef_plots, by = c("var" = ".variable", "k", "m")) |> 
    name_coef_indices(m, trial_col, stan_data) |> 
    select(var, covar, k, prior, likelihood, diagnosis, plot_obj) 
}

get_coef_sd_powerscale_table_data <- function(coef_ps_sense, prior_crcr_coef_sd, crcr_coef_sd, prior_crcr_pfs_coef_sd, crcr_pfs_coef_sd, stan_data) {
  coef_plots <- bind_rows(
    crcr_covar_trial_sd = bind_rows(prior = prior_crcr_coef_sd, posterior = crcr_coef_sd, .id = "fit_type") |> rename(.exp_value = crcr_covar_trial_sd),
    covar_trial_sd = bind_rows(prior = prior_crcr_pfs_coef_sd, posterior = crcr_pfs_coef_sd, .id = "fit_type") |> rename(.exp_value = covar_trial_sd),
    .id = ".variable"
  ) |>
    nest(coef_sd_data = !c(.variable, m)) |> 
    mutate(plot_obj = map(coef_sd_data, \(d) plot_coef(d, .exp_value))) 
  
  coef_ps_sense |> 
    filter(str_detect(variable, "trial_sd")) |> 
    tidyr::separate_wider_regex(variable, c(var = ".+", r"{\[}", m = r"{\d+}", "]")) |> 
    mutate(m = as.integer(m)) |> 
    left_join(coef_plots, by = c("m", "var" = ".variable")) |> 
    name_coef_indices(m, NULL, stan_data) |> 
    select(var, covar, prior, likelihood, diagnosis, plot_obj) 
}

get_covar_trial_sd <- function(res, stan_data) {
  spread_rvars(res, covar_trial_sd[m]) |> 
    name_coef_indices(m, NULL, stan_data)
}

get_joint_gng_prob <- function(mpfs_res_data, pfs6_res_data, orr_res_data, mpfs_cutoffs, pfs6_cutoffs, orr_cutoffs) {
  cutoffs <- bind_rows(mpfs = mpfs_cutoffs, pfs6 = pfs6_cutoffs, orr = orr_cutoffs, .id = "endpoint")

  bind_rows(
    mpfs = select(mpfs_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_median_pfs) |> 
      mutate(endpoint_forecast_val = weeks_to_months(endpoint_forecast_val)),
    pfs6 = select(pfs6_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_pfs6),
    orr = select(orr_res_data, model_type, fit_type, trial, endpoint_forecast_val = forecast_trial_subpop_orr),
    .id = "endpoint"
  ) |> 
    left_join(cutoffs, by = "endpoint") |> 
    pivot_wider(id_cols = c(model_type, fit_type, trial), names_from = endpoint, values_from = c(lrv, tv, endpoint_forecast_val)) |> 
    mutate(
      p_tv = Pr(endpoint_forecast_val_mpfs > tv_mpfs & endpoint_forecast_val_pfs6 > tv_pfs6 & endpoint_forecast_val_orr > tv_orr), 
      p_lrv = Pr(endpoint_forecast_val_mpfs > lrv_mpfs & endpoint_forecast_val_pfs6 > lrv_pfs6 & endpoint_forecast_val_orr > lrv_orr)
    ) 
}

plot_coef <- function(d, xvar) {
  ggplot(d) + 
    stat_slab(aes(xdist = {{ xvar }}, color = fit_type), fill = NA, linewidth = 2, show.legend = TRUE) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 2) +
    scale_color_discrete("", label = str_to_title, type = AZ_palette, aesthetic = c("color", "fill")) +
    labs(x = "", y = "") +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 25), legend.text = element_text(size = 20)) +
    coord_cartesian(xlim = c(0, 4)) +
    theme(legend.position = "bottom") +
    NULL
}

get_coef_powerscale_table_data <- function(coef_ps_sense, prior_crcr_coef, crcr_coef, prior_crcr_pfs_coef, crcr_pfs_coef, stan_data, trial_col = NULL) {
  coef_plots <- bind_rows(
    bind_rows(prior = prior_crcr_coef, posterior = crcr_coef, .id = "fit_type"),
    bind_rows(prior = prior_crcr_pfs_coef, posterior = crcr_pfs_coef, .id = "fit_type")
  ) |>
    filter(fct_match(.variable, c("crcr_covar_effect", "covar_effect", "crcr_tumor_stim_pop_coef", "tumor_stim_pop_coef"))) |> 
    nest(coef_data = !c(.variable, m, k)) |> 
    mutate(plot_obj = map(coef_data, \(d) plot_coef(d, .exp_value))) 
  
  coef_ps_sense |> 
    filter(!str_detect(variable, "trial_sd")) |> 
    tidyr::separate_wider_regex(variable, c(var = ".+", r"{\[}", trial = r"{\d+}", ",", m = r"{\d+}", r"{,?}", k = r"{(?:\d+)?}", ".*")) |> 
    mutate(
      across(c(m, k, trial), as.integer),
      k = factor(k, levels = 1:2, labels = c("Non-response", "Response"))
    ) |> 
    left_join(coef_plots, by = c("var" = ".variable", "k", "m")) |> 
    name_coef_indices(m, trial_col, stan_data) |> 
    select(var, covar, k, prior, likelihood, diagnosis, plot_obj) 
}

get_coef_sd_powerscale_table_data <- function(coef_ps_sense, prior_crcr_coef_sd, crcr_coef_sd, prior_crcr_pfs_coef_sd, crcr_pfs_coef_sd, stan_data) {
  coef_plots <- bind_rows(
    crcr_covar_trial_sd = bind_rows(prior = prior_crcr_coef_sd, posterior = crcr_coef_sd, .id = "fit_type") |> rename(.exp_value = crcr_covar_trial_sd),
    covar_trial_sd = bind_rows(prior = prior_crcr_pfs_coef_sd, posterior = crcr_pfs_coef_sd, .id = "fit_type") |> rename(.exp_value = covar_trial_sd),
    .id = ".variable"
  ) |>
    nest(coef_sd_data = !c(.variable, m)) |> 
    mutate(plot_obj = map(coef_sd_data, \(d) plot_coef(d, .exp_value))) 
  
  coef_ps_sense |> 
    filter(str_detect(variable, "trial_sd")) |> 
    tidyr::separate_wider_regex(variable, c(var = ".+", r"{\[}", m = r"{\d+}", "]")) |> 
    mutate(m = as.integer(m)) |> 
    left_join(coef_plots, by = c("m", "var" = ".variable")) |> 
    name_coef_indices(m, NULL, stan_data) |> 
    select(var, covar, prior, likelihood, diagnosis, plot_obj) 
}

# Tumor analysis plots #####

plot_prior_post_dens <- function(res_data, param = .value, normalize = "all") {
  res_data |> 
    ggplot(aes(xdist = {{ param }}, color = fit_type)) +
    stat_slab(aes(fill = fit_type), alpha = 0.25, normalize = normalize) +
    stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single"), .width = c(0.5, 0.8, 0.99)) +
    # stat_spike(at = "median") +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title, aesthetics = c("fill", "color")) +
    scale_y_continuous("", breaks = NULL) +
    NULL
}

plot_prior_post_hist <- function(res_data, param = .value, normalize = "all", ...) {
  res_data |> 
    ggplot(aes(xdist = {{ param }}, color = fit_type)) +
    stat_histinterval(aes(fill = fit_type), alpha = 0.25, normalize = normalize, ...) +
    # stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single"), .width = c(0.5, 0.8, 0.99)) +
    # stat_spike(at = "median") +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title, aesthetics = c("fill", "color")) +
    scale_y_continuous("", breaks = NULL) +
    NULL
}

plot_corr_decay <- function(res_data, param = .value) {
  res_data |> 
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = {{ param }}, fill = fit_type, color = fit_type), alpha = 0.25, .width = c(0.5, 0.8), linewidth = 0.5) +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title, aesthetics = c("fill", "color")) +
    labs(x = "Week", y = "Correlation") +
    NULL
}

plot_dynamics <- function(data, var, expect_rvar = TRUE, na.rm = FALSE) {
  pobj <- ggplot(data, aes(week))
  
  if (expect_rvar) {
    pobj <- pobj + stat_lineribbon(aes(ydist = {{ var }}, fill = stage), na.rm = na.rm, alpha = 0.25, linewidth = 0.5, .width = c(0.5, 0.8))
  } else {
    pobj <- pobj + stat_lineribbon(aes(y = {{ var }}, fill = stage), na.rm = na.rm, alpha = 0.25, .width = c(0.5, 0.8))
  }
  
  pobj +
    geom_point(aes(y = mmsumdiam), color = AZ_navy, size = 1.5, alpha = 0.75) +
    scale_fill_discrete("Stage", type = AZ_palette, label = c("obs" = "Observed", "forecast" = "Forecast")) +
    facet_wrap(vars(i), scales = "free") +
    NULL
}

plot_level_rates <- function(res_data) {
  ggplot(res_data) +
    geom_lineribbon(aes(x, .value, ymin = .lower, ymax = .upper, color = fit_type, fill = fit_type, group = .width), alpha = 0.25, step = "hv") +
    scale_color_discrete("", type = AZ_palette, aesthetics = c("color", "fill"), label = str_to_title) +
    scale_x_continuous("") + 
    labs(y = "") +
    facet_wrap(vars(.variable), scales = "free") + #, labeller = labeller(.variable = \(l) str_remove(l, "log_"))) +
    # coord_cartesian(xlim = c(0, 10)) +
    NULL
}

plot_level_decrease_prop <- function(res_data) {
  ggplot(res_data) +
    geom_lineribbon(aes(x, .value_exp, ymin = .lower, ymax = .upper, color = fit_type, fill = fit_type, group = .width), alpha = 0.25, step = "hv") +
    scale_color_discrete("", type = AZ_palette, aesthetics = c("color", "fill"), label = str_to_title) +
    scale_x_continuous("", breaks = seq(-1, 1, 0.2)) +
    scale_y_continuous("", breaks = NULL) +
    NULL
}

plot_confusion_matrix <- function(data, recorded, calculated, p, n) {
  data |> 
    mutate(
      nvar = {{ n }},
      pvar = {{ p }},
      n_label = if (!is_null(nvar)) str_glue("(n={ nvar })") else "",
      size_label = str_glue("{round(pvar, 3)}
                             {n_label}")
    ) |>  
    ggplot(aes(x = {{ recorded }}, y = {{ calculated }})) +
    geom_tile(aes(fill = {{ p }}), alpha = 0.5, color = "white", linewidth = 0.5) +
    geom_text(aes(label = size_label), color = AZ_darkpurple, size = 3) +
    scale_fill_gradient(low = AZ_turquoise, high = AZ_pink, name = "Proportion") +
    scale_x_discrete(limits = fct_rev) +
    facet_wrap(vars(trial), labeller = labeller(.default = str_to_upper)) +
    coord_fixed() +
    theme_minimal() +
    theme(panel.grid.major = element_blank()) +
    NULL
}

plot_ssls_coef <- function(res_data, name_var = n) {
  res_data |> 
    ggplot(aes(y = {{ name_var }})) +
    stat_pointinterval(aes(xdist = .value, color = fit_type), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 0)
}

# Prepare data for plotting
prepare_recist_plot_data <- function(data) {
  data |>
    select(!matches("((forecast|rep)_)?recist")) |> 
    group_by(i) |> 
    mutate(succ_week = lead(week, default = max(week) + 1)) |> 
    rowwise() |> 
    reframe(across(everything()), week = seq(week, succ_week - 1)) |> 
    pivot_longer(c(CR, PR, SD, PD), names_to = "recist_cat", values_to = "prob") |>
    mutate(recist_cat = factor(recist_cat, levels = c("CR", "PR", "SD", "PD")))
}

plot_recist_predictions <- function(data, 
                                   x_breaks = months_to_weeks(seq(0, 48, 12)),
                                   y_label = "Posterior Probability",
                                   caption = "Bar height represents probability; colors show RECIST categories;\nColored points represent observed RECIST.") {
  # Create the plot
  data |>
    prepare_recist_plot_data() |>
    ggplot(aes(week, y = prob)) +
    geom_col(aes(alpha = stage, fill = recist_cat), position = "fill", width = 1.01, linewidth = 0) +
    geom_vline(aes(xintercept = week), linetype = "dashed", 
               data = \(d) filter(d, fct_match(stage, "obs")) |> group_by(i) |> slice_max(week)) +
    geom_vline(aes(xintercept = pfs), linetype = "dashed", color = "white", 
               data = \(d) filter(d, !right_censored)) +
    geom_point(aes(y = 0.9, fill = response, shape = "obs"), 
               color = "black", size = 2, 
               show.legend = c(fill = TRUE, color = FALSE),
               data = \(d) filter(d, fct_match(stage, "obs")) |> group_by(i, succ_week) |> slice_min(week)) +
    geom_point(aes(y = 0.8, fill = det_response, shape = "target"), 
               color = "black", size = 2.5, 
               show.legend = c(fill = TRUE, color = FALSE),
               data = \(d) filter(d, fct_match(stage, "obs")) |> group_by(i, succ_week) |> slice_min(week)) +
    scale_x_continuous("Months", breaks = x_breaks, label = label_weeks_to_months) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_manual(values = c("CR" = AZ_green, "PR" = AZ_turquoise, "SD" = AZ_gold, "PD" = AZ_plum), 
                      name = "RECIST Category", 
                      aesthetics = c("color", "fill")) +
    scale_alpha_manual("", values = c(obs = 0.25, forecast = 0.5), labels = c(obs = "Observed", forecast = "Forecast")) +
    scale_shape_manual("", values = c(obs = 21, target = 23), labels = c(obs = "Observed", target = "Target Lesions Only")) +
    labs(
      y = y_label,
      caption = caption
    ) + 
    facet_wrap(vars(i)) +
    theme(legend.position = "bottom") + 
    NULL
}

plot_pfs_ppc <- function(data) {
  data |> 
    ggplot(aes(pfs + interval_censored + 1)) +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_x_continuous("Recorded PFS [Months]", breaks = months_to_weeks(seq(0, 48, 6)), label = label_weeks_to_months) +
    scale_y_continuous("Posterior PFS [Months]", breaks = months_to_weeks(seq(0, 48, 6)), label = label_weeks_to_months) +
    stat_pointinterval(aes(ydist = spop_target_pfs, color = event_type), .width = 0.8, alpha = 0.5, linewidth = 1, size = 0.5) +
    # geom_label_repel(aes(y = median(spop_pfs), label = i), size = 2.5) +
    scale_color_discrete("Event Type", labels = c(death = "Death", target_pd = "Target PD", nontarget_pd = "Non-target PD"), type = AZ_palette) +
    labs(caption = "Restricted to uncensored patients.") +
    facet_wrap(vars(trial), scales = "free", labeller = labeller(trial = str_to_upper)) +
    theme(legend.position = "bottom") + 
    NULL
}
  
# Distogram #######

# First, create a helper function for the row-adding adjustment
# This ensures consistency between Stat and Geom implementations
add_extra_visualization_row <- function(data) {
  
}

bin_dist <- function(data, ..., breaks, add_right_boundary = TRUE) {
  data |> 
    reframe(
      x = breaks[-length(breaks)], 
      across(c(...), \(d) rvar_sample_hist(d, breaks))
    ) %>%
    bind_rows(if (add_right_boundary) filter(., x == nth(breaks, -2)) |> mutate(x = last(breaks)))
}

# StatDistogram implementation
StatDistogram <- ggproto(
  "StatDistogram", ggdist:::StatLineribbon,
  default_params = c(ggdist:::StatLineribbon$default_params, freq = TRUE),
  
  compute_panel = function(self, data, scales, orientation = "horizontal", ...) {
    # Call parent method
    result <- ggproto_parent(ggdist:::StatLineribbon, self)$compute_panel(
      data, scales, orientation = orientation, ...
    )
    
    return(result)
  },
  
  setup_params = function(self, data, params) {
    params <- ggproto_parent(ggdist:::StatLineribbon, self)$setup_params(data, params)
    
    if (is_empty(params$breaks)) {
      params$breaks <- breaks_fixed(data$x, width = 30)
    }
    
    if (is_null(params$freq) || inherits(params$freq, "waiver")) {
      params$freq <- TRUE
    }
    
    return(params)
  },
  
  setup_data = function(self, data, params) {
    data <- data |> 
      group_by(across(!ydist)) |> 
      reframe(
        x = params$breaks[-length(params$breaks)], 
        ydist = rvar_sample_hist(ydist, params$breaks, freq = params$freq)
      ) %>%
      bind_rows(filter(., x == nth(params$breaks, -2)) |> mutate(x = last(params$breaks)))
      
    
    # Add the extra row before calling parent's setup_data
    # data <- add_extra_visualization_row(data)
    
    # Call the parent's setup_data
    data <- ggproto_parent(ggdist:::StatLineribbon, self)$setup_data(data, params)
    
    return(data)
  } 
)

# stat_distogram function
stat_distogram <- function(mapping = NULL, data = NULL,
                           geom = "lineribbon",
                           position = "identity",
                           ...,
                           step = "hv",
                           breaks = waiver(),
                           .width = c(0.5, 0.8, 0.95),
                           point_interval = "median_qi",
                           orientation = NA,
                           na.rm = FALSE,
                           show.legend = NA,
                           inherit.aes = TRUE,
                           freq = waiver()) {
  layer(
    stat = StatDistogram,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      step = step,
      breaks = breaks,
      .width = .width,
      point_interval = point_interval,
      orientation = orientation,
      na.rm = na.rm,
      freq = freq,
      ...
    )
  )
}

# nolint end: object_usage_linter
