get_lfo_cutoff_days <- function(first_cutoff_date, last_date, first_cutoff_day_idx, days_increment) {
  len <- time_length(last_date - first_cutoff_date, unit = "days") %/% days_increment + 1
  n <- seq(len) 
  
  tibble(
    n,
    cutoff_date = accumulate(n[-len], \(prev, n) prev + days(days_increment), .init = first_cutoff_date),
    cutoff_calendar_day = first_cutoff_day_idx + time_length(cutoff_date - first_cutoff_date, unit = "days"),
  ) 
}

lfo <- function(stan_data, model, cutoffs, output_path, basename, output_timestamp = FALSE, refit_n = min(cutoffs$n), k_threshold = 0.7, lean = TRUE, verbose = FALSE, exact = FALSE) {
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
      iter_warmup = 300, iter_sampling = 500, parallel_chains = 4, adapt_delta = 0.9,
      init = create_crcr_pfs_initializer(.),
      output_dir = file.path(output_path, "fit"), output_basename = str_glue("{basename}-{refit_n}"),
      timestamp = output_timestamp
    ) |> 
    lfo_log_lik() |> 
    mutate(refit_n, n = n + refit_n - 1) |> 
    left_join(select(remaining_cutoffs, n, cutoff_date, cutoff_calendar_day), by = "n")
  
  if (lean) {
    psis_results <- psis_results |> 
      select(!c(psis_obj, lwt))
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
    return(bind_rows(psis_results, lfo(model, stan_data, cutoffs, output_path, basename, output_timestamp, refit_n = min(next_cutoffs$n), k_threshold, lean, verbose)))
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

lfo_pointwise_log_lik <- function(res) {
  res |> 
    spread_rvars(oos_log_lik[n, m], oos_pfs_log_lik[n, m], oos_crcr_log_lik[n, m]) %>% {
      inner_join(
        filter(., m == max(m)) |> select(!m), 
        filter(., n == 1) |> select(!n), 
        by = c("n" = "m"), suffix = c("", "_log_ratio"))
    } |> 
    mutate(
      fit = map(min_rank(n), \(nr) if (nr == 1) res),
      across(ends_with("log_ratio"), lag),
      across(ends_with("log_lik"), \(oos_ll) map_dbl(oos_ll, \(ll) log_mean_exp(posterior::draws_of(ll))), .names = "original_E_log_lik_{.col}")
    ) |> 
    rename(original_E_log_lik = original_E_log_lik_oos_log_lik) |> 
    rename_with(\(col) str_replace(col, "_log_lik_oos_(pfs|crcr)_log_lik", r"{_\1_log_lik}"), starts_with("original_E_log_lik"))
}

lfo_log_lik <- function(res) {
  res |> 
    lfo_pointwise_log_lik() |>
    mutate(
      across(
        ends_with("log_ratio"), 
        \(oos_lr) map_if(oos_lr, \(lr) !is.na(lr), \(lr) suppressWarnings(loo::psis(posterior::draws_of(lr))), .else = \(lr) NA),
        .names = "psis_obj_{.col}"
      ),
    ) |> 
    rename_with(\(col) str_replace(col, "_oos((?:_pfs|_crcr)?)_log_lik_log_ratio", r"{\1}"), starts_with("psis_obj")) |>
    mutate(
      across(starts_with("psis_obj"), \(po) unlist(map_if(po, \(p) !is_na(p), loo::pareto_k_values)), .names = "k_{.col}"),
      across(starts_with("psis_obj"), \(po) map_if(psis_obj, \(p) !is_na(p), \(o) weights(o, normalize = TRUE)[, 1]), .names = "lwt_{.col}")
    ) |> 
    rename_with(\(col) str_replace(col, "_psis_obj((?:_pfs|_crcr)?)", r"{\1}"), c(starts_with("lwt_"), starts_with("k_"))) |>
    mutate(
      approx_E_log_lik = map2_vec(oos_log_lik, lwt, \(lr, lw) log_sum_exp(lw + posterior::draws_of(lr))),
      approx_E_log_lik_pfs = map2_vec(oos_pfs_log_lik, lwt_pfs, \(lr, lw) log_sum_exp(lw + posterior::draws_of(lr))),
      approx_E_log_lik_crcr = map2_vec(oos_crcr_log_lik, lwt_pfs, \(lr, lw) log_sum_exp(lw + posterior::draws_of(lr))),
    )
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
