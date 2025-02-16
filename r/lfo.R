get_lfo_cutoff_days <- function(first_cutoff_date, last_date, first_cutoff_day_idx, days_increment) {
  len <- time_length(last_date - first_cutoff_date, unit = "days") %/% days_increment + 1
  n <- seq(len) 
  
  tibble(
    n,
    cutoff_date = accumulate(n[-len], \(prev, n) prev + days(days_increment), .init = first_cutoff_date),
    cutoff_calendar_day = first_cutoff_day_idx + time_length(cutoff_date - first_cutoff_date, unit = "days"),
  ) 
}

lfo <- function(model, stan_data, cutoffs, output_path, basename, output_timestamp = FALSE, refit_n = min(cutoffs$n), k_threshold = 0.7, lean = TRUE, verbose = FALSE, exact = FALSE) {
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
    return(bind_rows(psis_results, lfo(model, stan_data, cutoffs, output_path, basename,  output_timestamp, refit_n = min(next_cutoffs$n), k_threshold, lean, verbose)))
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
    spread_rvars(oos_log_lik[n]) |>
    rowwise() |> 
    mutate(E_log_lik = log_mean_exp(posterior::draws_of(oos_log_lik[1]))) |> 
    ungroup()
}

lfo_log_lik <- function(res) {
  res |> 
    lfo_pointwise_log_lik() |>
    mutate(
      log_ratio = map(oos_log_lik, posterior::rvar_sum), 
      psis_obj = map(log_ratio, \(lr) suppressWarnings(loo::psis(posterior::draws_of(lr)))),
      k = lag(map_dbl(psis_obj, loo::pareto_k_values), default = NA_real_),
      lwt = lag(map(psis_obj, \(o) weights(o, normalize = TRUE)[, 1]), default = NA),
      approx_E_log_lik = map2_dbl(oos_log_lik, lwt, \(ll, lw) if (is_null(lw)) NA_real_ else log_sum_exp(lw + posterior::draws_of(ll)))
    )  
}
