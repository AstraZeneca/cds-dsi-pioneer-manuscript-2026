get_mpfs_gng_table <- function(res_data, model_type_names) {
  res_data |>
    mutate(
      model_type = fct_recode(
        model_type,
        !!!set_names(names(model_type_names), model_type_names)
      )
    ) |>
    gt(rowname_col = "model_type") |>
    tab_header(md("**Forecast Median PFS**")) |>
    # cols_move_to_start(c(.lower, forecast_trial_median_pfs, .upper, p_lrv, p_tv)) |>
    fmt_number(!model_type, decimals = 2) |>
    fmt_percent(c(p_lrv, p_tv)) |>
    tab_spanner("80% Credible Interval", c(.lower, post_median, .upper)) |>
    tab_spanner("G/NG Probability", starts_with("p_")) |>
    cols_align("center") |>
    cols_label(
      .lower = "Lower Bound",
      post_median = "Median",
      .upper = "Upper Bound",
      p_lrv = "LRV", #md(r"{$\Pr[\textrm{mPFS} > \textrm{LRV}]$}"),
      p_tv = "TV", #md(r"{$\Pr[\textrm{mPFS} > \textrm{TV}]$}"),
    )
}

get_pfs6_gng_table <- function(res_data, model_type_names) {
  res_data |>
    mutate(
      model_type = fct_recode(
        model_type,
        !!!set_names(names(model_type_names), model_type_names)
      )
    ) |>
    gt(rowname_col = "model_type") |>
    tab_header(md("**Forecast PFS6**")) |>
    fmt_number(!model_type, decimals = 2) |>
    fmt_percent(contains("_p_")) |>
    tab_spanner("80% Credible Interval", c(.lower, post_median, .upper)) |>
    tab_spanner("G/NG Probability", contains("_p_")) |>
    cols_align("center") |>
    cols_label(
      .lower = "Lower Bound",
      post_median = "Median",
      .upper = "Upper Bound",
      p_lrv = "LRV",
      p_tv = "TV",
    )
}

get_orr_gng_table <- function(res_data, model_type_names) {
  res_data |>
    mutate(
      model_type = fct_recode(
        model_type,
        !!!set_names(names(model_type_names), model_type_names)
      )
    ) |>
    gt(rowname_col = "model_type") |>
    tab_header(md("**Forecast ORR**")) |>
    fmt_number(!model_type, decimals = 2) |>
    fmt_percent(starts_with("p_")) |>
    tab_spanner("80% Credible Interval", c(.lower, post_median, .upper)) |>
    tab_spanner("G/NG Probability", starts_with("p_")) |>
    cols_align("center") |>
    cols_label(
      .lower = "Lower Bound",
      post_median = "Median",
      .upper = "Upper Bound",
      p_lrv = "LRV",
      p_tv = "TV",
    )
}

get_mpfs_table_data <- function(res_data, lrv, tv) {
  res_data |>
    filter(fct_match(fit_type, "posterior"), fct_match(trial, "lung")) |>
    mutate(
      across(ends_with("median_pfs"), weeks_to_months),
      p_lrv = Pr(forecast_trial_median_pfs > lrv),
      p_tv = Pr(forecast_trial_median_pfs > tv)
    ) |>
    point_interval(forecast_trial_median_pfs, .width = 0.8) |>
    select(
      model_type,
      .lower,
      post_median = forecast_trial_median_pfs,
      .upper,
      p_lrv,
      p_tv
    )
}

get_pfs6_table_data <- function(res_data, lrv, tv) {
  res_data |>
    filter(fct_match(fit_type, "posterior"), fct_match(trial, "lung")) |>
    mutate(
      # across(starts_with("forecast"), lst(p_lrv = \(x) Pr(x > lrv), p_tv = \(x) Pr(x > tv)))
      p_lrv = Pr(forecast_trial_pfs6 > lrv),
      p_tv = Pr(forecast_trial_pfs6 > tv),
    ) |>
    point_interval(forecast_trial_pfs6, .width = 0.8) |>
    select(
      model_type,
      .lower,
      post_median = forecast_trial_pfs6,
      .upper,
      p_lrv,
      p_tv
    )
}

get_orr_table_data <- function(res_data, lrv, tv) {
  res_data |>
    filter(fct_match(fit_type, "posterior"), fct_match(trial, "lung")) |>
    mutate(
      p_lrv = Pr(forecast_trial_subpop_orr > lrv),
      p_tv = Pr(forecast_trial_subpop_orr > tv)
    ) |>
    point_interval(forecast_trial_orr, .width = 0.8) |>
    select(
      model_type,
      .lower,
      post_median = forecast_trial_orr,
      .upper,
      p_lrv,
      p_tv
    )
}

get_combined_results_table <- function(res_data, model_type_names) {
  res_data |>
    select(!post_median) |>
    mutate(
      endpoint = fct_recode(
        endpoint,
        "mPFS" = "mpfs",
        "PFS6" = "pfs6",
        "ORR" = "orr"
      )
    ) |>
    gt(rowname_col = "endpoint", groupname_col = "model_type") |>
    row_group_order(c(
      "separate",
      "multilevel",
      "full_multilevel",
      "stacked"
    )) |>
    text_transform(
      \(x) {
        fct_recode(
          unlist(x),
          !!!set_names(names(model_type_names), model_type_names)
        )
      },
      locations = cells_row_groups()
    ) |>
    text_replace("^0%", "<1%", locations = cells_body(c(p_lrv, p_tv))) |>
    cols_merge(
      c(.lower, .upper),
      pattern = "[{1}, {2}]",
      rows = !fct_match(endpoint, "mPFS")
    ) |>
    cols_merge(
      c(.lower, .upper),
      pattern = "[{1}, {2}] months",
      rows = fct_match(endpoint, "mPFS")
    ) |>
    fmt_number() |>
    fmt_percent(c(p_lrv, p_tv), decimals = 0) |>
    tab_footnote(
      "All intervals are the 80% credible intervals.",
      cells_column_labels(.lower)
    ) |>
    tab_footnote(
      "Probability > Least Reference Value. mPFS = 6 months, PFS6 = 50%, ORR = 30%.",
      cells_column_labels(p_lrv)
    ) |>
    tab_footnote(
      "Probability > Target Value. mPFS = 9 months, PFS6 = 59%, ORR = 45%.",
      cells_column_labels(p_tv)
    ) |>
    cols_align("center") |>
    cols_label(
      .lower = "Forecasting",
      .upper = "Upper Bound",
      p_lrv = "LRV",
      p_tv = "TV"
    )
}

#' Create ORR summary table by estimate type and cohort
#'
#' @param estimate_type Either "sample" or "spop"
#' @param cohort Either "all", "first_line", or "part_e"
#' @param orr_data ORR data frame from all_tumor_ssls_orr_rvar_ctdna_jan26
#' @param cond_orr_data Conditional ORR data frame from all_tumor_ssls_cond_orr_rvar_ctdna_jan26
#' @return A gt table object with ORR summary statistics
create_orr_table <- function(
  estimate_type = c("sample", "spop"),
  cohort = c("all", "first_line", "part_e"),
  orr_data,
  cond_orr_data
) {
  estimate_type <- match.arg(estimate_type)
  cohort <- match.arg(cohort)
  col_prefix <- if (estimate_type == "sample") "sample" else "spop"

  cohort_var <- case_when(
    cohort == "all" ~ "all",
    cohort == "first_line" ~ "pdl1_naive",
    cohort == "part_e" ~ "part_e_pdl1"
  )

  cohort_label <- case_when(
    cohort == "all" ~ "All Patients",
    cohort == "first_line" ~ "First Line Patients",
    cohort == "part_e" ~ "Part E"
  )

  bind_rows(
    orr_data |> rename(orr = !!str_c(col_prefix, "_target_orr")),
    cond_orr_data |> rename(orr = !!str_c("cond_", col_prefix, "_target_orr"))
  ) |>
    filter(fct_match(dco, "jan26")) |>
    prepare_pdl1_and_trial_info() |>
    filter(fct_match(variable, cohort_var), fct_match(trial, "sclc"), fct_match(fit_type, "posterior")) |>
    point_interval(orr, .width = c(0.9)) |>
    select(cond_group_name, orr, .lower, .upper) |>
    gt(rowname_col = NULL) |>
    tab_header(
      title = str_c("ORR Posterior Summary - ", cohort_label),
    ) |>
    cols_label(
      orr = "Median ORR",
      .lower = "5%",
      .upper = "95%",
      cond_group_name = ""
    ) |>
    tab_spanner(
      label = "Percentiles",
      columns = c(".lower", ".upper")
    ) |>
    fmt_percent(
      columns = c(orr, .lower, .upper),
      decimals = 2
    )
}

#' Create median survival summary table by estimate type and cohort
#'
#' @param estimate_type Either "sample" or "spop"
#' @param cohort Either "all", "first_line", or "part_e"
#' @param quant_data Quantile data frame (PFS or OS)
#' @param cond_quant_data Conditional quantile data frame (PFS or OS)
#' @param endpoint Either "pfs" or "os"
#' @return A gt table object with median survival summary statistics
create_median_survival_table <- function(
  estimate_type = c("sample", "spop"),
  cohort = c("all", "first_line", "part_e"),
  quant_data,
  cond_quant_data,
  endpoint = c("pfs", "os")
) {
  estimate_type <- match.arg(estimate_type)
  cohort <- match.arg(cohort)
  endpoint <- match.arg(endpoint)
  col_prefix <- if (estimate_type == "sample") "sample" else "spop"

  if (endpoint == "pfs") {
    col_name <- str_c(col_prefix, "_pfs_quant")
    cond_col_name <- str_c("cond_", col_prefix, "_pfs_quant")
  } else {
    col_name <- str_c(col_prefix, "_os_quant")
    cond_col_name <- str_c("cond_", col_prefix, "_os_quant")
  }

  label <- if (endpoint == "pfs") "PFS" else "OS"

  cohort_var <- case_when(
    cohort == "all" ~ "all",
    cohort == "first_line" ~ "pdl1_naive",
    cohort == "part_e" ~ "part_e_pdl1"
  )

  cohort_label <- case_when(
    cohort == "all" ~ "All Patients",
    cohort == "first_line" ~ "First Line Patients",
    cohort == "part_e" ~ "Part E"
  )

  bind_rows(
    quant_data |> rename(median_val = !!col_name),
    cond_quant_data |> rename(median_val = !!cond_col_name)
  ) |>
    filter(quantile == 0.5, fct_match(dco, "jan26")) |>
    prepare_pdl1_and_trial_info() |>
    filter(fct_match(variable, cohort_var), fct_match(trial, "sclc"), fct_match(fit_type, "posterior")) |>
    point_interval(median_val, .width = c(0.9)) |>
    mutate(
      median_val = weeks_to_months(median_val),
      .lower = weeks_to_months(.lower),
      .upper = weeks_to_months(.upper)
    ) |>
    select(cond_group_name, median_val, .lower, .upper) |>
    gt(rowname_col = NULL) |>
    tab_header(
      title = str_c("Median ", label, " Posterior Summary - ", cohort_label),
    ) |>
    cols_label(
      median_val = str_c("Median ", label, " (months)"),
      .lower = "5%",
      .upper = "95%",
      cond_group_name = ""
    ) |>
    tab_spanner(
      label = "Percentiles",
      columns = c(".lower", ".upper")
    ) |>
    fmt_number(
      columns = c(median_val, .lower, .upper),
      decimals = 2
    )
}

#' @rdname create_median_survival_table
create_median_pfs_table <- function(estimate_type, cohort, pfs_quant_data, cond_pfs_quant_data) {
  create_median_survival_table(estimate_type, cohort, pfs_quant_data, cond_pfs_quant_data, endpoint = "pfs")
}

#' Create survival-n summary table by estimate type and cohort
#'
#' @param estimate_type Either "sample" or "spop"
#' @param cohort Either "all", "first_line", or "part_e"
#' @param surv_n_data Survival-n data frame (PFS or OS)
#' @param cond_surv_n_data Conditional survival-n data frame (PFS or OS)
#' @param endpoint Either "pfs" or "os"
#' @return A gt table object with survival-n summary statistics
create_survival_n_table <- function(
  estimate_type = c("sample", "spop"),
  cohort = c("all", "first_line", "part_e"),
  surv_n_data,
  cond_surv_n_data,
  endpoint = c("pfs", "os")
) {
  estimate_type <- match.arg(estimate_type)
  cohort <- match.arg(cohort)
  endpoint <- match.arg(endpoint)
  col_prefix <- if (estimate_type == "sample") "sample" else "spop"

  if (endpoint == "pfs") {
    col_name <- str_c(col_prefix, "_pfs_n")
    cond_col_name <- str_c("cond_", col_prefix, "_pfs_n")
  } else {
    col_name <- str_c(col_prefix, "_os_n")
    cond_col_name <- str_c("cond_", col_prefix, "_os_n")
  }

  label <- if (endpoint == "pfs") "PFS-n" else "OS-n"

  cohort_var <- case_when(
    cohort == "all" ~ "all",
    cohort == "first_line" ~ "pdl1_naive",
    cohort == "part_e" ~ "part_e_pdl1"
  )

  cohort_label <- case_when(
    cohort == "all" ~ "All Patients",
    cohort == "first_line" ~ "First Line Patients",
    cohort == "part_e" ~ "Part E"
  )

  bind_rows(
    surv_n_data |> rename(surv_n = !!col_name),
    cond_surv_n_data |> rename(surv_n = !!cond_col_name)
  ) |>
    filter(fct_match(dco, "jan26")) |>
    prepare_pdl1_and_trial_info() |>
    filter(fct_match(variable, cohort_var), fct_match(trial, "sclc"), fct_match(fit_type, "posterior"), timepoint %in% c(6, 9, 12, 18)) |>
    point_interval(surv_n, .width = c(0.9)) |>
    mutate(
      timepoint_label = str_c(timepoint, " Months")
    ) |>
    select(timepoint_label, cond_group_name, surv_n, .lower, .upper) |>
    group_by(timepoint_label) |>
    gt() |>
    tab_header(
      title = str_c(label, " Posterior Summary - ", cohort_label),
    ) |>
    cols_label(
      surv_n = str_c("Median ", label),
      .lower = "5%",
      .upper = "95%",
      cond_group_name = ""
    ) |>
    tab_style(
      style = cell_text(align = "center"),
      locations = cells_row_groups()
    ) |>
    tab_spanner(
      label = "Percentiles",
      columns = c(".lower", ".upper")
    ) |>
    fmt_percent(columns = c("surv_n", ".lower", ".upper"), decimals = 2)
}

#' @rdname create_survival_n_table
create_pfs_n_table <- function(estimate_type, cohort, pfs_n_data, cond_pfs_n_data) {
  create_survival_n_table(estimate_type, cohort, pfs_n_data, cond_pfs_n_data, endpoint = "pfs")
}
