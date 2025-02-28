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
     output_timestamp = FALSE, refit_n = min(cutoffs$n), k_threshold = 0.7, lean = TRUE, verbose = FALSE, exact = FALSE, iter_warmup = 300, iter_sampling = 500, ...) {
  if (verbose) {
    cat("Startin on:\n")
    print(cutoffs)
    cat("\n")
  }
  
  remaining_cutoffs <- cutoffs |> filter(n >= refit_n) 
  
  psis_results <- stan_data |>
    list_assign(cutoff_calendar_day = remaining_cutoffs$cutoff_calendar_day, n_cutoffs = nrow(remaining_cutoffs)) %>%
    sample_and_save(
      model,
      .,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling, parallel_chains = 4, adapt_delta = 0.9,
      init = create_crcr_pfs_initializer(.),
      output_dir = file.path(output_path, "fit"), output_basename = str_glue("{basename}-{refit_n}"),
      timestamp = output_timestamp, 
      ...
    ) |> 
    lfo_log_lik() |> 
    mutate(refit_n, n = n + refit_n - 1) |> 
    left_join(select(remaining_cutoffs, n, cutoff_date, cutoff_calendar_day), by = "n")
  
  if (lean) {
    psis_results <- psis_results |>
      select(!c(starts_with("psis"), starts_with("lwt"), fit))
  }
  
  next_cutoffs <- psis_results |> 
    filter(!map_lgl(k, is_null), map_dbl(k, max) > k_threshold | exact, n > refit_n) %>%
    semi_join(remaining_cutoffs, ., by = "n")
  
  if (verbose) {
    cat("LFO results:\n")
    print(psis_results)
    cat("\n")
  }
  
  if (nrow(next_cutoffs) > 0) {
    return(bind_rows(
      psis_results, 
      lfo(stan_data, model, cutoffs, output_path, basename, output_timestamp, refit_n = min(next_cutoffs$n), k_threshold, lean, verbose, exact, iter_warmup, iter_sampling, ...)
    ))
  } else {
    return(psis_results)
  }
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

lfo_log_lik <- function(res) {
  res |> 
    spread_rvars(patient_log_lik[n, m, i], patient_pfs_log_lik[n, m, i], patient_crcr_log_lik[n, m, i]) |>
    filter(map_lgl(patient_log_lik, \(l) any(l != 0))) %>% {
      inner_join(
        filter(., m == max(m)) |> select(!m), 
        filter(., n == 1) |> select(!n), 
        by = c("n" = "m", "i"), suffix = c("", "_log_ratio"))
    } |> 
    group_by(n) |> 
    summarize(across(matches("^patient(_pfs|_crcr)?_log_lik"), \(l) list(draws_of(l)))) |> 
    mutate(
      fit = map(min_rank(n), \(nr) if (nr == 1) res),
      across(matches("^patient(_pfs|_crcr)?_log_lik$"), \(l) map(l, \(ln) plyr::aaply(ln, 2, log_mean_exp)), .names = "mean_{.col}"),
      across(matches("^patient(_pfs|_crcr)?_log_lik_log_ratio$"), \(l) map(l, \(ln) suppressWarnings(loo::psis(ln))), .names = "psis_{.col}"), 
      across(
        starts_with("psis"), 
        lst(k = \(po) map(po, loo::pareto_k_values), lwt = \(po) map(po, \(pon) weights(pon, normalize = TRUE))), 
        .names = "{.fn}_{.col}"
      ),
      across(matches("^(psis|lwt|k)"), lag),
      dplyover::across2(
        matches("^patient(_pfs|_crcr)?_log_lik$"), matches("^lwt(_pfs|crcr)?"), 
        \(l, w) map2(l, w, \(ln, wn) if (!is_null(wn)) plyr::aaply(wn + ln, 2, log_mean_exp)), .names = "approx_mean_{xcol}"),
      across(matches("^(approx_)?mean"), \(m) map_dbl(m, sum), .names = "E_{.col}")
    ) |> 
    rename_with(\(n) str_replace_all(
      n, 
      c(r"{log_lik_log_ratio}" = "log_ratio",
        r"{(k|lwt)_psis_patient(_pfs|_crcr)?_log_ratio}" = r"{\1\2}", 
        r"{^psis_patient(_pfs|_crcr)?_log_ratio}" = r"{psis\1}",
        r"{E_(approx_)?mean}" = r"{\1E}")
    ))  
}

clean_lfo_results <- function(lfo_res) {
  lfo_res |> 
    group_by(n) %>% 
    filter(if (has_name(., "refit_n")) min_rank(refit_n) == n() else TRUE) |> 
    ungroup() |> 
    mutate(
      E_log_lik = if_else(is.na(k), original_E_log_lik, approx_E_log_lik),
      E_log_lik_pfs = if_else(is.na(k_pfs), original_E_pfs_log_lik, approx_E_log_lik_pfs),
      E_log_lik_crcr = if_else(is.na(k_crcr), original_E_crcr_log_lik, approx_E_log_lik_crcr),
    )
}

lfo_stacking_weights <- function(...) {
  model_log_lik <- rlang::dots_list(..., .named = TRUE)

  model_log_lik |> 
    map_dfr(clean_lfo_results, .id = "model") |> 
    select(model, n, E_log_lik) |> 
    pivot_wider(names_from = model, values_from = E_log_lik) |> 
    select(!n) |> 
    as.matrix() |> 
    loo::stacking_weights() |> 
    c() |> 
    set_names(names(model_log_lik))
}
