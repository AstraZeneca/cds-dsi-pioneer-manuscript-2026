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

plot_crcr_baseline_hazard <- function(res_data, analysis_data) {
  plot_baseline_hazard(res_data, crcr_trial_lambda, k) +
    geom_rug(
      aes(x = confirmed_response_week),
      alpha = 0.5,
      data = analysis_data |>
        filter(!confirmed_response_censored) |>
        mutate(k = if_else(confirmed_response, "Response", "Non-response"))
    ) +
    facet_grid(vars(trial), vars(k), scales = "free_y")
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
  ggplot(res_data, aes(time, estimate)) +
    stat_lineribbon(aes(fill = fit_type), linewidth = 0, alpha = 0.25, .width = c(0.5, 0.8)) +
    geom_step(aes(y = estimate, linetype = "Observed"), direction = "vh", data = obs_cif_data) +
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
  select(res_data, !.value) |> 
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
  ggplot(res_data, aes(y = covar)) +
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