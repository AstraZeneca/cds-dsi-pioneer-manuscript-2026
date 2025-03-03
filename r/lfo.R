get_lfo_cutoff_days <- function(first_cutoff_date, last_date, first_cutoff_day_idx, days_increment) {
  len <- time_length(last_date - first_cutoff_date, unit = "days") %/% days_increment + 1
  n <- seq(len) 
  
  tibble(
    n,
    cutoff_date = accumulate(n[-len], \(prev, n) prev + days(days_increment), .init = first_cutoff_date),
    cutoff_calendar_day = first_cutoff_day_idx + time_length(cutoff_date - first_cutoff_date, unit = "days"),
  ) 
}

lfo <- function(
    stan_data, model, cutoffs, output_path, basename, 
    output_timestamp = FALSE, refit_n = min(cutoffs$n), 
    k_threshold = 0.7, lean = FALSE, verbose = FALSE, exact = FALSE, fit_only = FALSE, iter_warmup = 300, iter_sampling = 500, ...) {
  if (verbose) {
    cat("Startin on:\n")
    print(cutoffs)
    cat("\n")
  }
  
  remaining_cutoffs <- cutoffs |> filter(n >= refit_n) 
  
  fit <- stan_data |>
    list_assign(cutoff_calendar_day = remaining_cutoffs$cutoff_calendar_day, n_cutoffs = nrow(remaining_cutoffs)) %>%
    sample_and_save(
      model,
      .,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling, parallel_chains = 4, adapt_delta = 0.9,
      init = create_crcr_pfs_initializer(.),
      output_dir = file.path(output_path, "fit"), output_basename = str_glue("{basename}-{refit_n}"),
      timestamp = output_timestamp, 
      ...
    ) 
  
  if (fit_only) {
    return(fit)
  }
  
  psis_results <- fit |> 
    lfo_log_lik() |> 
    mutate(refit_n, n = n + refit_n - 1) |> 
    left_join(select(remaining_cutoffs, n, cutoff_date, cutoff_calendar_day), by = "n")
  
  if (lean) {
    psis_results <- psis_results |>
      select(n, contains("E_"))
  }
  
  next_cutoffs <- psis_results |> 
    filter(!is.na(k), k > k_threshold | exact, n > refit_n) %>%
    semi_join(remaining_cutoffs, ., by = "n")
  
  if (verbose) {
    cat("LFO results:\n")
    print(psis_results)
    cat("\n")
  }
  
  if (nrow(next_cutoffs) > 0) {
    return(bind_rows(
      psis_results, 
      lfo(
        stan_data, model, cutoffs, output_path, basename, output_timestamp, refit_n = min(next_cutoffs$n), 
        k_threshold, lean, verbose, exact, fit_only, iter_warmup, iter_sampling, ...
      )
    ))
  } else {
    return(psis_results)
  }
}

lfo_drop_bad_approx <- function(lfo_res) {
  lfo_res |> 
    group_by(n) %>% 
    filter(if (has_name(., "refit_n")) min_rank(refit_n) == n() else TRUE) |> 
    ungroup()  
}

# more stable than log(sum(exp(x))) 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

lfo_log_lik <- function(res, max_n = Inf, future_window = 1) {
  res |> 
    spread_rvars(patient_log_lik[n, m, i], patient_pfs_log_lik[n, m, i], patient_crcr_log_lik[n, m, i]) |>
    lfo_log_lik_rvar(max_n, future_window)
}

psis_resample <- function(l, w, recalc_full = FALSE) { #, negative_only = TRUE) {
  map2(l, w, function(ln, wn) { 
    if (!is_null(wn)) {
      plyr::aaply(ln, 2, \(lni) log_sum_exp(lni + wn * all(wn < 0))) 
    } else if (recalc_full) { 
      plyr::aaply(ln, 2, log_mean_exp)
    }
  })
}

lfo_log_lik_rvar <- function(log_lik_rvar, max_n = Inf, future_window = 1) {
  log_lik_rvar |>   
    # filter(map_lgl(patient_log_lik, \(l) any(l != 0)), m >= n) |>  
    filter(m >= n) |>
    group_by(n, m) |> 
    summarize(across(matches("^patient(_pfs|_crcr)?_log_lik"), \(l) list(draws_of(l))), .groups = "drop") |>
    (function(d) {
      inner_join(
        filter(d, m == max(m)) |> select(!m), # From n to max(m)
        filter(d, n == 1) |> select(!n),      # From 1 to n
        by = c("n" = "m"), suffix = c("", "_log_ratio")
      ) |> 
        left_join(
          mutate(d, m = m - future_window + 1) |> filter(n == m), 
          by = "n", suffix = c("", "_w")
        )
    })() |>
    filter(n <= max_n) |> 
    mutate(
      # fit = map(min_rank(n), \(nr) if (nr == 1) res),
      across(
        matches("^patient(_pfs|_crcr)?_log_lik(_w)?$"), 
        \(l) map(l, \(ln) plyr::aaply(ln[, colSums(ln) < 0], 2, \(lni) log_mean_exp(lni))), 
        .names = "mean_{.col}"
      ),
      across(
        matches("^patient(_pfs|_crcr)?_log_lik_log_ratio$"),
        \(l) map(l, \(ln) suppressWarnings(loo::psis(rowSums(ln)))), 
        .names = "psis_{.col}"
      ), 
      across(
        starts_with("psis"), 
        lst(k = \(po) map_dbl(po, loo::pareto_k_values), lwt = \(po) map(po, \(pon) weights(pon, normalize = TRUE)[, 1])), 
        .names = "{.fn}_{.col}"
      ),
      across(matches("^(psis|lwt|k)"), lag),
      
    ) |> 
    rename_with(\(n) str_replace_all(
      n, 
      c(r"{log_lik_log_ratio}" = "log_ratio",
        r"{(k|lwt)_psis_patient(_pfs|_crcr)?_log_ratio}" = r"{\1\2}", 
        r"{^psis_patient(_pfs|_crcr)?_log_ratio}" = r"{psis\1}")
    )) |>   
    mutate(
      dplyover::across2(matches("^patient(_pfs|_crcr)?_log_lik$"), matches("^lwt(_pfs|crcr)?"), psis_resample, .names = "approx_mean_{xcol}"),
      dplyover::across2(matches("^patient(_pfs|_crcr)?_log_lik_w$"), matches("^lwt(_pfs|crcr)?"), psis_resample, .names = "approx_mean_{xcol}"),
      across(matches("^(approx_)?mean"), \(m) map_dbl(m, \(mn) if (!is_null(mn)) sum(mn) else NA_real_), .names = "E_{.col}")
    ) |> 
    rename_with(\(n) str_replace(n, r"{E_(approx_)?mean}", r"{\1E}"))
}

redo_lfo_results <- function(lfo_res, lean = FALSE) {
  new_res <- lfo_res |> 
    lfo_drop_bad_approx() |> 
    arrange(n) |> 
    group_by(refit_n) |> 
    reframe(lfo_log_lik(first(fit), max_n = n())) |>   
    mutate(n = n + refit_n - 1)
  
  if (lean) {
    new_res <- new_res |>
      select(n, contains("E_"), k)
  }
  
  return(new_res)
}

clean_lfo_results <- function(lfo_res) {
  lfo_res |> 
    lfo_drop_bad_approx() |> 
    mutate(
      E_log_lik = if_else(is.na(k), E_patient_log_lik, approx_E_patient_log_lik),
      E_pfs_log_lik = if_else(is.na(k), E_patient_pfs_log_lik, approx_E_patient_pfs_log_lik),
      E_crcr_log_lik = if_else(is.na(k), E_patient_crcr_log_lik, approx_E_patient_crcr_log_lik),
      E_log_lik_w = if_else(is.na(k), E_patient_log_lik_w, approx_E_patient_log_lik_w),
      E_pfs_log_lik_w = if_else(is.na(k), E_patient_pfs_log_lik_w, approx_E_patient_pfs_log_lik_w),
      E_crcr_log_lik_w = if_else(is.na(k), E_patient_crcr_log_lik_w, approx_E_patient_crcr_log_lik_w),
    )
}

lfo_bootstrap_elpd <- function(lfo_res, n_bootstrap = 1000) {
  get_patient_subset_col <- function(ln, i) {
    n_early_patients <- length(i) - ncol(ln)
    ln[, discard(i, \(x) x < n_early_patients) - n_early_patients]
  } 
  
  origin_res <- lfo_res |> 
    lfo_drop_bad_approx() |> 
    arrange(n)
  
  n_lfo_patients <- ncol(first(origin_res$patient_log_lik))
  
  map_dfr(seq(n_bootstrap), function(b) {
    bootstrap_i <- sample(n_lfo_patients, n_lfo_patients, replace = TRUE)
   
    origin_res |>
      mutate(
        across(matches("^patient(_pfs|_crcr)?_log_lik(_w)?$"), \(l)  map(l, \(ln) get_patient_subset_col(ln, bootstrap_i))),
        dplyover::across2(
          matches("^patient(_pfs|_crcr)?_log_lik$"), matches("^lwt(_pfs|crcr)?"), \(l, lw) psis_resample(l, lw, recalc_full = TRUE), .names = "mean_{xcol}"
        ),
        dplyover::across2(
          matches("^patient(_pfs|_crcr)?_log_lik_w$"), matches("^lwt(_pfs|crcr)?"), \(l, lw) psis_resample(l, lw, recalc_full = TRUE), .names = "mean_{xcol}"
        ),
      ) |> 
      transmute(across(matches("^mean"), \(m) map_dbl(m, \(mn) if (!is_null(mn)) sum(mn) else NA_real_), .names = "E_{.col}")) |> 
      rename_with(\(n) str_replace(n, r"{E_(approx_)?mean_patient}", r"{\1E}")) |> 
      summarize(across(everything(), sum))
  }) |> 
    summarize(across(everything(), sd))
}

lfo_stacking_weights <- function(..., log_lik_var = E_log_lik) {
  model_log_lik <- rlang::dots_list(..., .named = TRUE)

  model_log_lik |> 
    map_dfr(clean_lfo_results, .id = "model") |> 
    select(model, n, {{ log_lik_var }}) |> 
    pivot_wider(names_from = model, values_from = {{ log_lik_var }}) |> 
    select(!n) |> 
    as.matrix() |> 
    loo::stacking_weights() |> 
    c() |> 
    set_names(names(model_log_lik))
}
