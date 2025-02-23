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

plot_cif <- function(res_data, obs_cif_data, time = time, estimate = estimate) {
  ggplot(res_data, aes({{ time }}, {{ estimate }})) +
    stat_lineribbon(aes(fill = fit_type), linewidth = 0, alpha = 0.25, .width = c(0.5, 0.8)) +
    geom_step(aes(x = time, y = estimate, linetype = "Observed"), direction = "vh", data = \(d) semi_join(obs_cif_data, d, by = "trial")) +
    scale_linetype_manual("", values = c(Observed = "dashed")) +
    labs(y = "CIF") +
    NULL
}

plot_cif2 <- function(res_data, obs_cif_data) {
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

plot_km <- function(res_data, obs_km_data, km_est, group = fit_type) {
  ggplot(res_data) +
    geom_step(aes(x = t, y = s, group = btype, color = btype), linewidth = 0.5, alpha = 0.5, data = \(d) semi_join(obs_km_data, d, by = "trial")) +
    stat_lineribbon(aes(x = t - 1, ydist = {{ km_est }}, fill = {{ group }}, alpha = {{ group }}), linewidth = 0, .width = 0.8) +
    labs(y = "Survival Probability") +
    guides(alpha = "none") + 
    theme(legend.position = "bottom")
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

plot_patient_timelines <- function(analysis_data) {
  analysis_data |> 
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
    geom_vline(xintercept = c(lubridate::ymd("2024-06-30"), lubridate::ymd("2024-10-31")), linetype = "dashed", color = AZ_platinum) +
    annotate("text", x = lubridate::ymd("2024-6-30") - days(40), y = 2, label = "DCO 1") +
    annotate("text", x = lubridate::ymd("2024-10-31") - days(40), y = 2, label = "DCO 2") +
    labs(x = "Calendar Time", y = "Patients") +
    scale_color_discrete("", label = c("FALSE" = "Progression", "TRUE" = "Censored"), type = AZ_palette) +
    scale_shape_manual(
      "", values = c("visit" = 124, "treat" = 5, "last" = 19), labels = c("visit" = "Visit", "treat" = "Treatment Start", "last" = "Last Visit")
    ) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "inside", legend.position.inside = c(0.25, 0.8))
}