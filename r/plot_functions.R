# nolint start: object_usage_linter

# Trial name labeller for plots
trial_labeller <- function(x) {
  case_match(
    x,
    "sclc" ~ "SCLC-01",
    .default = str_to_upper(x)
  )
}

prepare_pdl1_and_trial_info <- function(res_data) {
  res_data |>
    filter(
      (!fct_match(variable, "first_liners") | !fct_match(cond_group_name, "no")) &
      (!fct_match(variable, "parts") | fct_match(cond_group_name, "part_e"))
    ) |>
    mutate(
      trial = if ("trial" %in% names(res_data)) coalesce(trial, "sclc") else "sclc",
      across(c(cond_group_name, variable), \(l) coalesce(l, "all")),
      variable = fct_collapse(variable,
        "all" = c("all", "pdl1"),
        "pdl1_naive" = c("pdl1_naive", "first_liners"),
        "part_e_pdl1" = c("part_e_pdl1", "parts")
      ),
      cond_group_name = fct_collapse(cond_group_name, "All" = c("all", "yes", "part_e"), "PDL1 Low" = "low", "PDL1 High" = "hi")
    )
}

#' Plot outcome (e.g., ORR or Median PFS) by PDL1 status and trial for SCLC trial
#'
#' @param data Data frame with columns: outcome, cond_group_name, variable, trial, fit_type, etc.
#' @param outcome Unquoted column name for the outcome to plot (e.g., orr, median_pfs)
#' @param ... Additional arguments passed to ggplot2::stat_pointinterval
#' @return A ggplot object
plot_outcome_by_pdl1_and_trial <- function(res_data, outcome, .width = c(0.5, 0.9), ...) {
  res_data |>
    prepare_pdl1_and_trial_info() |>
    filter(fct_match(variable, c("all", "pdl1_naive", "part_e_pdl1")), fct_match(trial, "sclc")) |>
    ggplot() +
    stat_pointinterval(aes(xdist = {{ outcome }}, y = cond_group_name, color = fit_type), position = "dodge", .width = .width, ...) +
    scale_color_discrete("", type = AZ_palette, label = str_to_title) +
    facet_grid(
      vars(variable),
      vars(trial),
      scales = "free",
      space = "free",
      labeller = labeller(trial = trial_labeller, variable = c("all" = "All", "pdl1_naive" = "First Line", "part_e_pdl1" = "Part E"))
    ) +
    labs(caption = "Points show posterior median; inner bars show 50% credible intervals,\nouter bars show 90% credible intervals.") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_baseline_hazard <- function(res_data, lambda_var, ...) {
  ggplot(res_data, aes(t)) +
    geom_line(
      aes(y = {{ lambda_var }}, color = fit_type, group = str_c(fit_type, .draw)),
      alpha = 0.25,
      data = \(d) {
        unnest_rvars(d) |>
          ungroup() %>%
          semi_join(
            distinct(., fit_type, ..., .draw) |>
              group_by(fit_type, ...) |>
              slice_sample(n = 25),
            by = join_by(fit_type, ..., .draw)
          )
      }
    ) +
    stat_lineribbon(aes(ydist = {{ lambda_var }}, fill = fit_type, color = fit_type), alpha = 0.25, linewidth = 0, .width = 0.8) +
    labs(
      y = "Baseline Hazard",
      caption = "Ribbons represent 80% credible intervals; faint lines show 25 posterior draws."
    ) +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    )
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
      linewidth = 1,
      step = "hv",
      .width = c(0.5, 0.8)
    ) +
    geom_step(
      aes(y = count, linetype = "Observed"),
      direction = "vh",
      linewidth = 0.5,
      data = \(d) {
        analysis_data |>
          group_by(trial) |>
          reframe(t = conf_resp_hb[-length(conf_resp_hb)], count = sample_hist(confirmed_response_week, conf_resp_hb))
      }
    ) +
    scale_alpha_manual("", values = c(prior = 0.125, posterior = 0.25)) +
    scale_linetype_manual("", values = c(Observed = "dotted")) +
    labs(
      y = "Count",
      caption = "Ribbons represent 50% and 80% credible intervals; lines show posterior median."
    ) +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    ) +
    guides(alpha = "none") +
    NULL
}

plot_cif <- function(res_data, obs_cif_data) {
  ggplot(res_data) +
    stat_lineribbon(aes(x = t, ydist = trial_cif, fill = fit_type), linewidth = 0, alpha = 0.25, .width = c(0.5, 0.8)) +
    geom_step(aes(x = time, y = estimate, linetype = "Observed"), direction = "vh", data = \(d) semi_join(obs_cif_data, d, by = "trial")) +
    scale_linetype_manual("", values = c(Observed = "dashed")) +
    labs(
      y = "CIF",
      caption = "Ribbons represent 50% and 80% credible intervals around cumulative incidence."
    ) +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_hazard_ratio <- function(res_data) {
  ggplot(res_data, aes(t)) +
    stat_lineribbon(
      aes(ydist = bindist, color = fit_type, fill = fit_type, alpha = fit_type),
      linewidth = 1,
      step = "hv",
      .width = c(0.5, 0.8)
    ) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    scale_alpha_manual("", values = c(prior = 0.125, posterior = 0.25)) +
    labs(
      x = "Hazard Ratio",
      y = "Density",
      caption = "Ribbons represent 50% and 80% credible intervals; lines show posterior median."
    ) +
    guides(alpha = "none") +
    theme(
      legend.position = "bottom",
      strip.text.y.left = element_text(angle = 0),
      strip.placement = "outside",
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    ) +
    NULL
}

plot_crcr_covar_coef <- function(res_data) {
  res_data |>
    filter(fct_match(.variable, "crcr_covar_trial_coef")) |>
    select(!.value) |>
    pivot_wider(names_from = k, values_from = .rs_value) |>
    mutate(tumor = str_detect(covar, "tumor sizes"), hazard_ratio = Response / `Non-response`) |>
    ggplot(aes(y = covar)) +
    stat_pointinterval(aes(xdist = hazard_ratio, color = fit_type), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    labs(x = "Ratio of Competing Risk Hazards", y = "Parameter", caption = "Showing the posterior median, 50% CI, and 80% CI.") +
    coord_cartesian(clip = "off") +
    theme(legend.position = "bottom", strip.text.y = element_blank()) +
    NULL
}

plot_pfs_covar_coef <- function(res_data) {
  res_data |>
    filter(fct_match(.variable, "covar_trial_coef")) |>
    ggplot(aes(y = covar)) +
    stat_pointinterval(aes(xdist = .rs_value, color = fit_type), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    labs(x = "Exponential of Parameter", y = "Parameter", caption = "Showing the posterior median, 50% CI, and 80% CI.") +
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
      aes(yend = usubjid, x = {{ surv_interval_col }}, xend = {{ surv_interval_col }} + {{ ic_col }}),
      linewidth = 0.25,
      data = \(d) filter(d, !{{ rc_col }})
    ) +
    geom_segment(
      aes(yend = usubjid, x = {{ surv_interval_col }}, xend = max_week),
      linewidth = 0.25,
      linetype = "dashed",
      data = \(d) filter(d, {{ rc_col }})
    ) +
    geom_point(aes(x = {{ surv_interval_col }}), size = 0.5) +
    geom_point(aes(x = {{ surv_interval_col }} + {{ ic_col }}, color = {{ exit_color }}), size = 0.5, data = \(d) {
      filter(d, !{{ rc_col }})
    }) +
    scale_color_discrete("", type = AZ_palette, label = c("FALSE" = "Non-response", "TRUE" = "Response"), aesthetic = c("color", "fill")) +
    # scale_color_ramp_discrete() +
    labs(
      y = "",
      caption = "Horizontal bars represent 50% and 80% credible intervals for predicted survival times."
    ) +
    facet_grid(
      vars(trial),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = labeller(trial = c("endometrial" = "Endometrial", "lung" = "Lung"))
    ) +
    theme(
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "top",
      strip.text.y.left = element_text(angle = 0),
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    ) +
    guides(
      color_ramp = "none", # Remove the interval legend
      color = guide_legend("") # Keep only the confirmed_response legend
    ) +
    NULL
}

base_plot_km <- function(
  res_data,
  obs_km_data,
  km_est,
  analysis_data = NULL,
  group = fit_type,
  alpha_group = NULL,
  color_group = NULL,
  linewidth = 0,
  obs_line_color = NULL,
  endpoint = c("pfs", "os"),
  ...
) {
  endpoint <- match.arg(endpoint)

  # Handle defaults - check if NULL and assign
  group_quo <- enquo(group)
  alpha_group_quo <- enquo(alpha_group)
  color_group_quo <- enquo(color_group)

  if (rlang::quo_is_null(alpha_group_quo)) {
    alpha_group_quo <- group_quo
  }
  if (rlang::quo_is_null(color_group_quo)) {
    color_group_quo <- group_quo
  }

  pobj <- ggplot(res_data, aes(x = t - 1)) +
    stat_lineribbon(
      aes(
        dist = {{ km_est }},
        fill = !!group_quo,
        alpha = !!alpha_group_quo,
        group = !!group_quo,
      ),
      linewidth = linewidth,
      .width = 0.8,
      ...
    ) +
    # stat_lineribbon(
    #   aes(
    # dist = {{ km_est }},
    #   color = !!color_group_quo
    #   ),
    #   .width = 0, # only median line
    #   linewidth = linewidth,
    #   alpha = 1,
    #   fill = NA,
    #   show.legend = FALSE
    # ) +
    labs(
      y = "Survival Probability",
      caption = "Ribbons represent 80% credible intervals; lines show posterior median survival curves."
    ) +
    guides(alpha = "none") +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    )

  if (!is_null(obs_km_data)) {
    if (is.null(obs_line_color)) {
      # Default: color obs KM line by btype (upper/lower bound) — used by SCLC
      pobj <- pobj +
        geom_step(aes(x = t, y = s, group = btype, color = btype), linewidth = 0.75, alpha = 0.5, data = \(d) {
          semi_join(obs_km_data, d, by = "trial")
        }) +
        geom_point(
          aes(x = t, y = s, color = btype),
          shape = 3, size = 2, stroke = 0.8,
          data = \(d) semi_join(obs_km_data, d, by = "trial") |> filter(c > 0)
        )
    } else {
      # Fixed color obs KM line (used by pioneer — no btype distinction needed)
      pobj <- pobj +
        geom_step(
          aes(x = t, y = s),
          color = obs_line_color, linewidth = 0.75, linetype = "dashed",
          data = \(d) semi_join(obs_km_data, d, by = "trial")
        )
    }

    if (!is_null(analysis_data)) {
      has_censor_reason <- "censor_reason" %in% names(analysis_data)
      if (endpoint == "pfs") {
        if (has_censor_reason) {
          pobj <- pobj +
            geom_point(
              aes(x = t, y = s, color = censor_reason, shape = censor_reason),
              size = 2.5, alpha = 0.9,
              data = \(d) {
                semi_join(obs_km_data, d, by = "trial") |>
                  inner_join(
                    analysis_data |> filter(right_censored) |> select(pfs, censor_reason),
                    by = c("t" = "pfs"), relationship = "many-to-many"
                  )
              }
            ) +
            geom_point(
              aes(x = t, y = s), color = "grey30", shape = "o",
              size = 2, alpha = 0.7,
              data = \(d) {
                semi_join(obs_km_data, d, by = "trial") |>
                  inner_join(
                    analysis_data |> filter(!right_censored, !progression_before_death) |> select(pfs),
                    by = c("t" = "pfs"), relationship = "many-to-many"
                  )
              }
            ) +
            scale_shape_manual(
              "Censoring reason",
              values = c("Admin censored" = "|", "Dropout" = "+")
            )
        } else {
          pobj <- pobj +
            geom_point(aes(x = t, y = s, color = btype, shape = "censored"), size = 2, alpha = 0.7, data = \(d) {
              semi_join(obs_km_data, d, by = "trial") |>
                inner_join(analysis_data |> filter(right_censored) |> select(pfs), by = c("t" = "pfs"), relationship = "many-to-many")
            }) +
            geom_point(aes(x = t, y = s, color = btype, shape = "death"), size = 2, alpha = 0.7, data = \(d) {
              semi_join(obs_km_data, d, by = "trial") |>
                inner_join(
                  analysis_data |> filter(!right_censored, !progression_before_death) |> select(pfs),
                  by = c("t" = "pfs"),
                  relationship = "many-to-many"
                )
            }) +
            scale_shape_manual("", values = c(censored = "|", death = "o"), labels = c(censored = "Right Censored", death = "Death before PD"))
        }
      } else if (endpoint == "os") {
        os_tick_data <- \(d) {
          semi_join(obs_km_data, d, by = "trial") |>
            inner_join(analysis_data |> filter(os_censored) |> select(os_time), by = c("t" = "os_time"), relationship = "many-to-many")
        }
        if (is.null(obs_line_color)) {
          pobj <- pobj +
            geom_point(aes(x = t, y = s, color = btype), shape = "|", size = 2, alpha = 0.7, data = os_tick_data)
        } else {
          pobj <- pobj +
            geom_point(aes(x = t, y = s), color = obs_line_color, shape = "|", size = 2, alpha = 0.7, data = os_tick_data)
        }
        pobj <- pobj +
          scale_shape_manual("", values = c(censored = "|"), labels = c(censored = "Censored (Alive)"))
      }
    }
  }

  return(pobj)
}

plot_gng <- function(res_data, outcome, lrv_tv, model_type_names) {
  res_data |>
    ggplot(aes(y = data_cut)) +
    stat_interval(
      aes(xdist = {{ outcome }}, color = impute_type, color_ramp = after_stat(level)),
      position = "dodge",
      .width = c(0.6, 0.8)
    ) +
    geom_vline(xintercept = lrv_tv, linetype = "dashed") +
    scale_x_continuous("", sec.axis = sec_axis(identity, breaks = lrv_tv, labels = c("LRV", "TV"))) +
    scale_y_discrete("", labels = str_to_title) +
    scale_color_discrete("Sample", type = AZ_palette, labels = \(l) str_replace(l, "_", " ") |> str_to_title()) +
    scale_color_ramp_discrete(name = "Credible Intervals") +
    facet_grid(vars(model_type), switch = "y", labeller = labeller(model_type = model_type_names)) +
    labs(caption = "Horizontal bars represent 60% and 80% credible intervals around the median.") +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    )
}

plot_simple_gng <- function(
  res_data,
  outcome,
  color_col,
  lrv_tv,
  model_type_names,
  y_col = model_type,
  outcome_desc = "",
  decision_prob = c(0.2, 0.9)
) {
  res_data |>
    ggplot(aes(y = {{ y_col }})) +
    stat_interval(
      aes(xdist = {{ outcome }}, color = {{ color_col }}, color_ramp = after_stat(level)),
      position = "dodge",
      .width = c(0.6, 0.8)
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
    labs(caption = "Horizontal bars represent 60% and 80% credible intervals around the median.") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

add_dco_to_plot <- function(plot, label_offset_x = -days(40), label_offset_y = 2) {
  plot +
    geom_vline(
      xintercept = c(lubridate::ymd("2024-04-22"), lubridate::ymd("2024-07-31"), lubridate::ymd("2024-10-31")),
      linetype = "dashed",
      color = AZ_platinum
    ) +
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
      "",
      values = c("visit" = 124, "treat" = 5, "last" = 19),
      labels = c("visit" = "Visit", "treat" = "Treatment Start", "last" = "Last Visit")
    ) +
    theme(
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.25, 0.8)
    )

  add_dco_to_plot(plot)
}

# Tumor analysis plots #####

plot_prior_post_dens <- function(res_data, param = .value, normalize = "all") {
  res_data |>
    ggplot(aes(xdist = {{ param }}, color = fit_type)) +
    stat_slab(aes(fill = fit_type), alpha = 0.25, normalize = normalize) +
    stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single"), .width = c(0.5, 0.8, 0.99)) +
    # stat_spike(at = "median") +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_colour_discrete("", type = AZ_palette, label = str_to_title) +
    scale_y_continuous("", breaks = NULL) +
    labs(caption = "Points show posterior median; intervals represent 50%, 80%, and 99%\ncredible intervals.") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_prior_post_hist <- function(res_data, param = .value, normalize = "all", ...) {
  res_data |>
    ggplot(aes(xdist = {{ param }}, color = fit_type)) +
    stat_histinterval(aes(fill = fit_type), alpha = 0.25, normalize = normalize, ...) +
    # stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single"), .width = c(0.5, 0.8, 0.99)) +
    # stat_spike(at = "median") +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_colour_discrete("", type = AZ_palette, label = str_to_title) +
    scale_y_continuous("", breaks = NULL) +
    NULL
}

plot_corr_decay <- function(res_data, param = .value) {
  res_data |>
    ggplot(aes(t)) +
    stat_lineribbon(aes(ydist = {{ param }}, fill = fit_type, color = fit_type), alpha = 0.25, .width = c(0.5, 0.8), linewidth = 0.5) +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_colour_discrete("", type = AZ_palette, label = str_to_title) +
    labs(
      x = "Week",
      y = "Correlation",
      caption = "Ribbons represent 50% and 80% credible intervals; lines show posterior median."
    ) +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
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
    labs(caption = "Ribbons represent 50% and 80% credible intervals;\nlines show posterior median tumor dynamics.") +
    facet_wrap(vars(i), scales = "free") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_level_param <- function(res_data, param = .value) {
  res_data |>
    filter(!is.na({{ param }})) |>
    ggplot() +
    geom_lineribbon(
      aes(x, {{ param }}, ymin = .lower, ymax = .upper, color = fit_type, fill = fit_type, group = .width),
      alpha = 0.25,
      step = "hv",
      linewidth = 0
    ) +
    scale_color_discrete("", type = AZ_palette, label = str_to_title) +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_x_continuous("") +
    labs(
      y = "",
      breaks = NULL,
      caption = "Ribbons represent credible intervals; lines show posterior median."
    ) +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_level_rates <- function(res_data) {
  ggplot(res_data) +
    geom_lineribbon(
      aes(x, .value, ymin = .lower, ymax = .upper, color = fit_type, fill = fit_type, group = .width),
      alpha = 0.25,
      step = "hv"
    ) +
    scale_color_discrete("", type = AZ_palette, label = str_to_title) +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_x_continuous("") +
    labs(
      y = "",
      caption = "Ribbons represent credible intervals; lines show posterior median."
    ) +
    facet_wrap(vars(.variable), scales = "free") + #, labeller = labeller(.variable = \(l) str_remove(l, "log_"))) +
    # coord_cartesian(xlim = c(0, 10)) +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_level_decrease_prop <- function(res_data) {
  ggplot(res_data) +
    geom_lineribbon(
      aes(x, .rs_value, ymin = .lower, ymax = .upper, color = fit_type, fill = fit_type, group = .width),
      alpha = 0.25,
      step = "hv"
    ) +
    scale_color_discrete("", type = AZ_palette, label = str_to_title) +
    scale_fill_discrete("", type = AZ_palette, label = str_to_title) +
    scale_x_continuous("", breaks = seq(-1, 1, 0.2)) +
    scale_y_continuous("", breaks = NULL) +
    labs(caption = "Ribbons represent credible intervals; lines show posterior median.") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))) +
    NULL
}

plot_confusion_matrix <- function(data, recorded, calculated, p, n = NULL) {
  nq <- enquo(n)

  data |>
    mutate(
      nvar = if (quo_is_null(nq)) "Unknown" else !!nq,
      pvar = {{ p }},
      n_label = if (!quo_is_null(nq)) str_glue("(n={ nvar })") else "",
      size_label = if_else(
        n_label == "",
        format(round(pvar, 2), nsmall = 2),
        str_glue("{format(round(pvar, 2), nsmall = 2)}\n{n_label}")
      )
    ) |>
    ggplot(aes(x = {{ recorded }}, y = {{ calculated }})) +
    geom_tile(aes(fill = {{ recorded }}, alpha = {{ p }}), color = "white", linewidth = 1.5) +
    geom_text(aes(label = size_label), color = AZ_darkpurple, size = 5, lineheight = 0.9, fontface = "bold") +
    scale_fill_recist(guide = "none") +
    scale_alpha_continuous(range = c(0, 0.9), guide = "none") +
    scale_x_discrete(limits = fct_rev, drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 15, 5, "pt")
    ) +
    NULL
}

# Reusable function for OOS confusion matrix plot
plot_oos_confusion_matrix <- function(confusion_matrix_data, add_caption = TRUE) {
  p <- confusion_matrix_data |>
    mutate(mp = median(mean_pred)) |>
    plot_confusion_matrix(response, pred_response, mp) +
    labs(
      x = "Observed RECIST Response",
      y = "Predicted RECIST Response"
    )

  if (add_caption) {
    p <- p + labs(caption = "Cell proportions show sensitivity (recall) for each RECIST category:\nP(Predicted response | Observed response). Diagonal elements indicate correct prediction rates.")
  }

  return(p)
}

# Compute specificity matrix from confusion matrix data
compute_specificity_matrix <- function(confusion_matrix_data) {
  # Get all unique response categories
  categories <- sort(unique(confusion_matrix_data$response))

  # Use count rvar if available, otherwise fall back to mean_pred
  if ("count" %in% names(confusion_matrix_data) && inherits(confusion_matrix_data$count, "rvar")) {
    # For each category X, compute P(Predicted = Y | Observed ≠ X)
    # using the count rvar to properly handle uncertainty
    specificity_data <- map_dfr(categories, function(cat) {
      # Get all observations where response != cat
      not_cat_data <- confusion_matrix_data |>
        filter(response != cat)

      # Total count for obs != cat
      total_not_cat <- rvar_sum(not_cat_data$count)

      # For each prediction category, sum counts and normalize
      not_cat_data |>
        group_by(pred_response) |>
        summarise(
          pred_count = rvar_sum(count),
          .groups = "drop"
        ) |>
        mutate(
          mean_pred = median(pred_count / total_not_cat),
          not_response = cat
        ) |>
        select(not_response, pred_response, mean_pred)
    })
  } else {
    # Fallback: use mean_pred and cell_size (less accurate)
    if ("mean_pred" %in% names(confusion_matrix_data) && inherits(confusion_matrix_data$mean_pred, "rvar")) {
      confusion_matrix_data <- confusion_matrix_data |>
        mutate(mean_pred = median(mean_pred))
    }

    # Compute prevalence of each observed category
    if ("cell_size" %in% names(confusion_matrix_data)) {
      obs_prevalence <- confusion_matrix_data |>
        group_by(response) |>
        summarise(prevalence = first(cell_size) / sum(unique(cell_size)), .groups = "drop")
    } else {
      obs_prevalence <- tibble(
        response = categories,
        prevalence = 1 / length(categories)
      )
    }

    specificity_data <- map_dfr(categories, function(cat) {
      not_cat_data <- confusion_matrix_data |>
        filter(response != cat) |>
        left_join(obs_prevalence, by = "response")

      total_not_cat_prevalence <- sum(filter(obs_prevalence, response != cat)$prevalence)

      not_cat_data |>
        mutate(conditional_prevalence = prevalence / total_not_cat_prevalence) |>
        group_by(pred_response) |>
        summarise(
          mean_pred = sum(mean_pred * conditional_prevalence),
          .groups = "drop"
        ) |>
        mutate(not_response = cat) |>
        select(not_response, pred_response, mean_pred)
    })
  }

  # Ensure we have all combinations
  all_combos <- expand_grid(
    not_response = categories,
    pred_response = categories
  )

  specificity_data <- all_combos |>
    left_join(specificity_data, by = c("not_response", "pred_response")) |>
    mutate(mean_pred = replace_na(mean_pred, 0))

  specificity_data
}

# Plot specificity matrix
plot_oos_specificity_matrix <- function(confusion_matrix_data) {
  specificity_data <- compute_specificity_matrix(confusion_matrix_data)

  specificity_data |>
    mutate(
      mp = mean_pred,
      # Create factor with explicit levels to ensure correct ordering
      not_response = factor(not_response, levels = c("CR", "PR", "SD", "PD"))
    ) |>
    plot_confusion_matrix(not_response, pred_response, mp) +
    labs(
      x = "Observed RECIST Response",
      y = "Predicted RECIST Response",
      caption = "Cell proportions show P(Predicted response | Observed ≠ category).\nColumn sums to 1.0. Higher off-diagonal values indicate better specificity."
    ) +
    scale_x_discrete(labels = ~paste0("NOT\n", .x))
}

# Reusable function for OOS confusion matrix Sankey diagram
plot_oos_confusion_sankey <- function(confusion_matrix_data) {
  # Load and prepare data
  data <- confusion_matrix_data |>
    mutate(
      mp = median(mean_pred),
      # Create a flag for correct predictions
      correct = response == pred_response,
      pct_label = if_else(correct, paste0(round(mp * 100, 1), "%"), NA_character_)
    ) |>
    # Ensure RECIST factor levels are in order
    mutate(
      response = factor(response, levels = c("CR", "PR", "SD", "PD")),
      pred_response = factor(pred_response, levels = c("CR", "PR", "SD", "PD"))
    )

  # Create base Sankey plot
  p <- ggplot(data,
         aes(y = mp, axis1 = response, axis2 = pred_response)) +
    geom_alluvium(aes(fill = response, alpha = correct), width = 1/12) +
    geom_stratum(width = 1/12, fill = "white", color = "grey30", linewidth = 0.5) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5)

  # Extract actual alluvium positions from the plot
  plot_build <- ggplot_build(p)
  alluvium_data <- plot_build$data[[1]]  # First layer is geom_alluvium

  # Filter for correct predictions (alpha = 1.0 after scale transformation)
  # and get positions at x=1 and x=2 to calculate center
  label_data <- alluvium_data |>
    filter(alpha == 1) |>  # Correct predictions have alpha = 1.0
    group_by(group, x) |>
    summarize(
      y_mid = mean((ymin + ymax) / 2),
      fill = first(fill),
      .groups = "drop"
    ) |>
    pivot_wider(names_from = x, values_from = y_mid, names_prefix = "y_x") |>
    mutate(
      y_center = 0.25 * y_x1 + 0.75 * y_x2,  # 3/4 of the way to the right
      x = 1.75  # 3/4 of the way between 1 and 2
    )

  # Get the percentage labels from original data
  correct_data <- data |>
    filter(correct) |>
    arrange(response) |>
    mutate(pct_label = paste0(round(mp * 100, 1), "%"))

  # Match by order (both should be in same order after filtering)
  label_data <- label_data |>
    arrange(desc(y_center)) |>  # Order by y position (top to bottom)
    mutate(pct_label = correct_data$pct_label)

  # Add labels to the plot
  p +
    geom_text(
      data = label_data,
      aes(x = x, y = y_center, label = pct_label),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold",
      color = "white"
    ) +
    scale_x_discrete(limits = c("Recorded\nResponse", "Predicted\nResponse"), expand = c(0.15, 0.05)) +
    scale_fill_recist() +
    scale_alpha_manual(
      values = c("TRUE" = 0.8, "FALSE" = 0.4),
      guide = "none"
    ) +
    labs(
      caption = "Flow width represents prediction proportions; percentages show correct prediction rates"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
}

plot_ssls_coef <- function(res_data, name_var = n, ...) {
  res_data |>
    ggplot(aes(y = {{ name_var }})) +
    stat_pointinterval(aes(xdist = .value, ...), point_size = 1, position = "dodge", .width = c(0.5, 0.8)) +
    geom_vline(xintercept = 0) +
    labs(caption = "Points show posterior median; inner bars show 50% credible intervals,\nouter bars show 80% credible intervals.") +
    theme(plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6)))
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

plot_recist_predictions <- function(
  data,
  x_breaks = months_to_weeks(seq(0, 48, 12)),
  y_label = "Posterior Probability",
  caption = "Bar height represents probability; colors show RECIST categories;\nColored points represent observed RECIST."
) {
  # Create the plot
  data |>
    prepare_recist_plot_data() |>
    ggplot(aes(week, y = prob)) +
    geom_col(aes(alpha = stage, fill = recist_cat), position = "fill", width = 1.01, linewidth = 0) +
    geom_vline(aes(xintercept = week), linetype = "dashed", data = \(d) {
      filter(d, fct_match(stage, "obs")) |> group_by(i) |> slice_max(week)
    }) +
    geom_vline(aes(xintercept = pfs), linetype = "dashed", color = "white", data = \(d) filter(d, !right_censored)) +
    geom_point(
      aes(y = 0.9, fill = response, shape = "obs"),
      color = "black",
      size = 2,
      show.legend = c(fill = TRUE, color = FALSE),
      data = \(d) filter(d, fct_match(stage, "obs")) |> group_by(i, succ_week) |> slice_min(week)
    ) +
    geom_point(
      aes(y = 0.8, fill = det_response, shape = "target"),
      color = "black",
      size = 2.5,
      show.legend = c(fill = TRUE, color = FALSE),
      data = \(d) filter(d, fct_match(stage, "obs")) |> group_by(i, succ_week) |> slice_min(week)
    ) +
    scale_x_continuous("Months", breaks = x_breaks, label = label_weeks_to_months) +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    scale_fill_recist() +
    scale_color_recist() +
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

plot_pfs_ppc <- function(data, pfs_var = spop_target_pfs, label_patients = FALSE) {
  plot_obj <- data |>
    ggplot(aes(pfs)) +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_x_continuous("Recorded PFS [Months]", breaks = months_to_weeks(seq(0, 48, 6)), label = label_weeks_to_months) +
    scale_y_continuous("Posterior PFS [Months]", breaks = months_to_weeks(seq(0, 48, 6)), label = label_weeks_to_months) +
    stat_pointinterval(aes(ydist = {{ pfs_var }}, color = event_type), .width = 0.8, alpha = 0.5, linewidth = 1, size = 0.5) +
    scale_color_discrete(
      "Event Type",
      labels = c(death = "Death", target_pd = "Target PD", nontarget_pd = "Non-target PD"),
      type = AZ_palette
    ) +
    labs(caption = "Points show posterior median with 80% credible intervals.\nRestricted to uncensored patients.") +
    facet_wrap(vars(trial), scales = "free", labeller = labeller(trial = trial_labeller)) +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
    ) +
    NULL

  if (label_patients) {
    plot_obj <- plot_obj +
      geom_label_repel(aes(y = median({{ pfs_var }}), label = i), size = 2.5) +
      NULL
  }

  return(plot_obj)
}

# Distogram #######

# First, create a helper function for the row-adding adjustment
# This ensures consistency between Stat and Geom implementations
add_extra_visualization_row <- function(data) {}

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
  "StatDistogram",
  ggdist:::StatLineribbon,
  default_params = c(ggdist:::StatLineribbon$default_params, freq = TRUE),

  compute_panel = function(self, data, scales, orientation = "horizontal", ...) {
    # Call parent method
    result <- ggproto_parent(ggdist:::StatLineribbon, self)$compute_panel(
      data,
      scales,
      orientation = orientation,
      ...
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
stat_distogram <- function(
  mapping = NULL,
  data = NULL,
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
  freq = waiver()
) {
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

#' Plot LFO ELPD difference comparing models against a baseline
#'
#' @param lfo_data Data frame with LFO results containing columns: model, n, m, E_log_lik_w (already filtered to desired models)
#' @param baseline_model Character string specifying which model to use as baseline
#' @param model_labels Named character vector for custom model labels (e.g., c(ctdna = "Covariates: All"))
#' @param add_vline Logical, whether to add a vertical line at x=0 (default: TRUE)
#' @param add_baseline_label Logical, whether to add a text label for the baseline model (default: TRUE)
#' @param baseline_label_text Character string for the baseline label text (optional, auto-derived from model_labels if available)
#' @param caption Character string for the plot caption
#' @return A ggplot object
plot_lfo_elpd_diff <- function(
  lfo_data,
  baseline_model,
  model_labels = NULL,
  add_vline = TRUE,
  add_baseline_label = TRUE,
  baseline_label_text = NULL,
  caption = "The cross bars indicate the standard errors and twice the standard errors of the mean elpd difference."
) {
  # Prepare the data
  plot_data <- lfo_data |>
    mutate(baseline_model = baseline_model) |>
    (\(.) left_join(., ., by = c("baseline_model" = "model", "n", "m"), suffix = c("", "_base"), relationship = "many-to-one"))() |>
    filter(model != baseline_model) |>
    mutate(E_log_lik_w_diff = E_log_lik_w - E_log_lik_w_base) |>
    group_by(model) |>
    summarize(
      across(c(E_log_lik_w_diff), lst(mean = \(ll) mean(ll, na.rm = TRUE), se = \(ll) sd(ll, na.rm = TRUE) / sqrt(n()))),
      .groups = "drop"
    ) |>
    pivot_longer(!model, names_to = c("metric", ".value"), names_pattern = "(.+)_(mean|se)$")

  # Calculate the maximum label width for dynamic margin adjustment
  if (!is.null(model_labels)) {
    label_texts <- model_labels[plot_data$model]
    label_texts <- label_texts[!is.na(label_texts)]
  } else {
    label_texts <- plot_data$model
  }

  # Calculate max characters and adjust left margin accordingly
  max_chars <- max(nchar(label_texts))
  # Estimate margin in lines: roughly 1 line per 10 characters, with minimum of 8
  left_margin_lines <- max(8, ceiling(max_chars / 10))

  # Calculate baseline label text and its width for right margin adjustment
  baseline_text_chars <- 0
  if (add_baseline_label) {
    if (is.null(baseline_label_text)) {
      # Try to get label from model_labels
      if (!is.null(model_labels) && !is.null(names(model_labels))) {
        # Look up the baseline model in the named vector
        temp_baseline_label <- model_labels[baseline_model]
        if (!is.na(temp_baseline_label)) {
          baseline_label_text <- paste("Baseline:", temp_baseline_label)
        } else {
          baseline_label_text <- paste("Baseline:\n", baseline_model)
        }
      } else {
        # Otherwise just use the baseline model name
        baseline_label_text <- paste("Baseline:", baseline_model)
      }
    }
    # Calculate the max line width in the baseline label text
    baseline_lines <- strsplit(baseline_label_text, "\n")[[1]]
    baseline_text_chars <- max(nchar(baseline_lines))
  }

  # Calculate right margin based on baseline text width
  # Estimate roughly 1 line per 8 characters for right margin (more generous), with minimum of 5.5
  right_margin_lines <- max(5.5, ceiling(baseline_text_chars / 8))

  # Start building the plot
  p <- ggplot(plot_data, aes(mean, model))

  # Add vertical line if requested
  if (add_vline) {
    p <- p + geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.8)
  }

  # Add baseline label if requested
  if (add_baseline_label) {
    p <- p + annotate("text", x = 0, y = Inf, label = baseline_label_text, hjust = 0.5, vjust = 1.5, size = 3.5, color = "gray30")
  }

  # Add crossbars
  p <- p +
    geom_crossbar(aes(xmin = mean - se, xmax = mean + se), fill = AZ_gold, alpha = 0.5, width = 0.25) +
    geom_crossbar(aes(xmin = mean - 2 * se, xmax = mean + 2 * se), fill = AZ_gold, alpha = 0.5, width = 0.25)

  # Calculate x-axis limits with extra space on the right for the baseline text
  x_min <- min(plot_data$mean - 2 * plot_data$se, na.rm = TRUE)
  x_max <- max(plot_data$mean + 2 * plot_data$se, na.rm = TRUE)

  # If adding a vertical line at zero, ensure zero is included in the axis
  if (add_vline) {
    x_min <- min(x_min, 0)
    x_max <- max(x_max, 0)
  }

  x_range <- x_max - x_min

  # Add extra space on the right: proportional to baseline text length
  # Use at least 20% of the range, or more if baseline text is long
  extra_right_space <- max(0.2 * x_range, baseline_text_chars * 0.01 * x_range)
  xlim <- c(x_min - 0.05 * x_range, x_max + extra_right_space)

  # Add scales
  if (!is.null(model_labels)) {
    p <- p + scale_y_discrete("", labels = model_labels)
  } else {
    p <- p + scale_y_discrete("")
  }

  p <- p +
    scale_x_continuous("Expected Log Predictive Density Difference") +
    labs(caption = caption) +
    coord_cartesian(xlim = xlim, clip = "off") +
    theme(
      plot.margin = margin(t = 5.5, r = right_margin_lines * 5.5, b = 5.5, l = left_margin_lines * 5.5, unit = "pt")
    )

  return(p)
}

# nolint end: object_usage_linter

plot_competing_risks_cif <- function(
  stan_data,
  draws_cif,
  cif_prefix = c("spop", "sample"),
  time_step = 4L,
  x_breaks_months = seq(0, 48, by = 6),
  trials = NULL
) {
  cif_prefix <- match.arg(cif_prefix)
  trial_id   <- as.integer(stan_data$patient_trial)
  n_trials   <- max(trial_id)
  trials     <- trials %||% seq_len(n_trials)
  trial_names <- trial_labeller(levels(stan_data$patient_trial) %||% as.character(seq_len(n_trials)))
  max_t      <- stan_data$max_all_t
  x_breaks   <- months_to_weeks(x_breaks_months)

  # ── Observed AJ CIF ──────────────────────────────────────────────────────
  obs <- tibble(
    trial   = trial_id,
    cens_01 = stan_data$ms_censored_01, t_01 = stan_data$ms_time_01,
    cens_02 = stan_data$ms_censored_02, t_02 = stan_data$ms_time_02,
    fs      = stan_data$ms_final_state,  t_03 = stan_data$ms_time_03
  ) |> mutate(
    cr_cause = case_when(
      fs == 3L ~ 3L, cens_01 == 0L ~ 1L, cens_02 == 0L ~ 2L, TRUE ~ 0L
    ),
    cr_time = case_when(
      cr_cause == 3L ~ t_03, cr_cause == 1L ~ t_01,
      cr_cause == 2L ~ t_02, TRUE ~ t_01
    )
  )

  cause_labels <- c("0\u21921 Progression", "0\u21922 On-trial death", "0\u21923 Dropout")

  aj_obs <- map_dfr(trials, function(tr) {
    sub <- filter(obs, trial == tr)
    cif <- cmprsk::cuminc(sub$cr_time, sub$cr_cause, cencode = 0)
    map_dfr(1:3, function(cause) {
      key <- paste0("1 ", cause)
      if (key %in% names(cif))
        tibble(time = cif[[key]]$time, est = cif[[key]]$est, cause = cause, trial = tr)
    })
  }) |> mutate(
    trial_name = trial_names[trial],
    cause_lbl  = factor(cause_labels[cause], levels = cause_labels)
  )

  # Observed state 0 retention (KM for any exit)
  obs_state0 <- map_dfr(trials, function(tr) {
    sub <- filter(obs, trial == tr)
    fit <- survival::survfit(survival::Surv(cr_time, cr_cause != 0) ~ 1, data = sub)
    tibble(time = fit$time, s0 = fit$surv, trial_name = trial_names[tr])
  })

  # ── Model CIF from Stan draws ────────────────────────────────────────────
  draws_mat <- posterior::as_draws_matrix(draws_cif)
  time_grid <- seq(time_step, max_t, by = time_step)

  model_cif <- map_dfr(trials, function(tr) {
    map_dfr(time_grid, function(t) {
      t_idx <- t + 1L
      extract_cif <- function(cause_id) {
        col <- paste0(cif_prefix, "_cif_0", cause_id, "[", tr, ",", t_idx, "]")
        if (!col %in% colnames(draws_mat)) return(tibble(med = NA_real_, lo = NA_real_, hi = NA_real_))
        x <- draws_mat[, col]
        q <- quantile(x, c(0.10, 0.50, 0.90))
        tibble(med = q[[2]], lo = q[[1]], hi = q[[3]])
      }
      bind_cols(
        tibble(t = t, trial = tr),
        rename_with(extract_cif(1), ~paste0("cif_01_", .)),
        rename_with(extract_cif(2), ~paste0("cif_02_", .)),
        rename_with(extract_cif(3), ~paste0("cif_03_", .))
      )
    })
  })

  model_cif_long <- model_cif |>
    pivot_longer(
      cols = starts_with("cif_"),
      names_to = c("cause", ".value"),
      names_pattern = "cif_(\\d+)_(\\w+)"
    ) |>
    mutate(
      cause      = as.integer(cause),
      trial_name = trial_names[trial],
      cause_lbl  = factor(cause_labels[cause], levels = cause_labels)
    )

  # Model state 0 retention
  model_state0 <- model_cif |>
    mutate(
      s0_med = 1 - cif_01_med - cif_02_med - cif_03_med,
      s0_lo  = pmax(0, 1 - cif_01_hi  - cif_02_hi  - cif_03_hi),
      s0_hi  = pmin(1, 1 - cif_01_lo  - cif_02_lo  - cif_03_lo),
      trial_name = trial_names[trial]
    )

  # ── Colors ────────────────────────────────────────────────────────────────
  cause_colors <- c(
    "0\u21921 Progression"    = AZ_navy,
    "0\u21922 On-trial death" = AZ_pink,
    "0\u21923 Dropout"        = AZ_turquoise
  )

  strip_theme <- theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = AZ_navy),
    strip.text = element_text(color = "white", face = "bold"),
    plot.caption = element_text(hjust = 0, size = rel(0.9), margin = margin(t = 6))
  )

  # ── Panel 1: CIF ─────────────────────────────────────────────────────────
  p_cif <- ggplot() +
    geom_ribbon(
      data = model_cif_long,
      aes(x = t, ymin = lo, ymax = hi, fill = cause_lbl),
      alpha = 0.2
    ) +
    geom_line(
      data = model_cif_long,
      aes(x = t, y = med, color = cause_lbl),
      linewidth = 0.8
    ) +
    geom_step(
      data = aj_obs,
      aes(x = time, y = est, color = cause_lbl),
      linetype = "dashed", linewidth = 0.9
    ) +
    facet_wrap(~trial_name, ncol = n_trials) +
    scale_color_manual(values = cause_colors, name = NULL) +
    scale_fill_manual(values = cause_colors, name = NULL) +
    scale_x_continuous(breaks = x_breaks, labels = label_weeks_to_months) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1), breaks = seq(0, 1, 0.2)
    ) +
    labs(
      title = "Competing risks CIF from state 0",
      x = "Months from baseline", y = "Cumulative incidence",
      caption = "Solid line/ribbon = posterior predictive median and 80% CI. Dashed = observed Aalen-Johansen."
    ) +
    theme_bw() +
    strip_theme

  # ── Panel 2: State 0 retention ───────────────────────────────────────────
  p_state0 <- ggplot() +
    geom_ribbon(
      data = model_state0,
      aes(x = t, ymin = s0_lo, ymax = s0_hi),
      fill = AZ_platinum, alpha = 0.4
    ) +
    geom_line(
      data = model_state0,
      aes(x = t, y = s0_med),
      color = AZ_navy, linewidth = 0.8
    ) +
    geom_step(
      data = obs_state0,
      aes(x = time, y = s0),
      color = "#F0AB00", linewidth = 0.9, linetype = "dashed"
    ) +
    facet_wrap(~trial_name, ncol = n_trials) +
    scale_x_continuous(breaks = x_breaks, labels = label_weeks_to_months) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, 1), breaks = seq(0, 1, 0.2)
    ) +
    labs(
      title = "Proportion still in state 0 (event-free)",
      x = "Months from baseline", y = "P(still in state 0)",
      caption = "Solid = model. Dashed = observed KM. Gap indicates the 0\u21921 hazard is too diffuse."
    ) +
    theme_bw() +
    strip_theme

  p_cif / p_state0 +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

#' Table of observed competing risks event counts from state 0
#'
#' @param stan_data Stan data list with patient_trial, ms_censored_01, etc.
#' @return A gt table summarising event counts by trial and cause.
