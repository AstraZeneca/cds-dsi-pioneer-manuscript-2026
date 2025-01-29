get_mpfs_gng_table <- function(res_data, model_type_names) {
  res_data |> 
    mutate(model_type = fct_recode(model_type, !!!set_names(names(model_type_names), model_type_names))) |> 
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
    mutate(model_type = fct_recode(model_type, !!!set_names(names(model_type_names), model_type_names))) |> 
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
    mutate(model_type = fct_recode(model_type, !!!set_names(names(model_type_names), model_type_names))) |> 
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
    select(model_type, .lower, post_median = forecast_trial_median_pfs, .upper, p_lrv, p_tv)
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
    select(model_type, .lower, post_median = forecast_trial_pfs6, .upper, p_lrv, p_tv)
}

get_orr_table_data <- function(res_data, lrv, tv) {
  res_data |> 
    filter(fct_match(fit_type, "posterior"), fct_match(trial, "lung")) |> 
    mutate(p_lrv = Pr(forecast_trial_subpop_orr > lrv), p_tv = Pr(forecast_trial_subpop_orr > tv)) |> 
    point_interval(forecast_trial_orr, .width = 0.8) |> 
    select(model_type, .lower, post_median = forecast_trial_orr, .upper, p_lrv, p_tv)
}

get_combined_results_table <- function(res_data, model_type_names) {
  res_data |> 
    select(!post_median) |> 
    mutate(endpoint = fct_recode(endpoint, "mPFS" = "mpfs", "PFS6" = "pfs6", "ORR" = "orr")) |> 
    gt(rowname_col = "endpoint", groupname_col = "model_type") |> 
    row_group_order(c("separate", "multilevel", "full_multilevel", "stacked")) |>
    text_transform(\(x) fct_recode(unlist(x), !!!set_names(names(model_type_names), model_type_names)), locations = cells_row_groups()) |> 
    text_replace("^0%", "<1%", locations = cells_body(c(p_lrv, p_tv))) |> 
    cols_merge(c(.lower, .upper), pattern = "[{1}, {2}]", rows = !fct_match(endpoint, "mPFS")) |> 
    cols_merge(c(.lower, .upper), pattern = "[{1}, {2}] months", rows = fct_match(endpoint, "mPFS")) |> 
    fmt_number() |> 
    fmt_percent(c(p_lrv, p_tv), decimals = 0) |>
    tab_footnote("All intervals are the 80% credible intervals.", cells_column_labels(.lower)) |> 
    tab_footnote("Probability > Least Reference Value. mPFS = 6 months, PFS6 = 50%, ORR = 30%.", cells_column_labels(p_lrv)) |> 
    tab_footnote("Probability > Target Value. mPFS = 9 months, PFS6 = 59%, ORR = 45%.", cells_column_labels(p_tv)) |> 
    cols_align("center") |> 
    cols_label(
      .lower = "Forecasting",
      .upper = "Upper Bound",
      p_lrv = "LRV",
      p_tv = "TV"
    )
}