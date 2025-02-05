get_lfo_cutoff_days <- function(first_cutoff_date, last_date, first_cutoff_day_idx, days_increment) {
  len <- time_length(last_date - first_cutoff_date, unit = "days") %/% days_increment + 1
  n <- seq(len) 
  
  tibble(
    n,
    cutoff_date = accumulate(n[-len], \(prev, n) prev + days(days_increment), .init = first_cutoff_date),
    cutoff_calendar_day = first_cutoff_day_idx + time_length(cutoff_date - first_cutoff_date, unit = "days"),
  ) 
}

lfo <- function(model, stan_data, cutoffs, basename, k_threshold = 0.7, lean = TRUE) {
  refit_n <- min(cutoffs$n)
  
  psis_results <- stan_data |>
    list_assign(cutoff_calendar_day = cutoffs$cutoff_calendar_day, n_oos_log_lik = nrow(cutoffs)) %>%
    sample_and_save(
      model,
      .,
      iter_warmup = 300, iter_sampling = 500, parallel_chains = 4, threads_per_chain = 4, adapt_delta = 0.9,
      init = create_crcr_pfs_initializer(.),
      output_dir = file.path(output_path, "fit"), output_basename = str_glue("{basename}-{refit_n}"),
      timestamp = fit_output_timestamp
    ) |> 
    lfo_log_lik() |> 
    mutate(refit_n, n = n + refit_n - 1) |> 
    left_join(select(cutoffs, n, cutoff_date, cutoff_calendar_day), by = "n")
  
  if (lean) {
    psis_results <- psis_results |> 
    select(n, refit_n, contains("log_lik"), k)
  }
  
  next_cutoffs <- psis_results |> 
    filter(!is.na(k), k > k_threshold, n > refit_n) %>%
    semi_join(cutoffs, ., by = "n")
  
  if (nrow(next_cutoffs) > 0) {
    return(bind_rows(psis_results, lfo(model, stan_data, next_cutoffs, basename, k_threshold, ...)))
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
