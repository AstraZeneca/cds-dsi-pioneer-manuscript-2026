

lfo <- function(model, N, L, M, k_threshold) {
  refits <- NULL 
  ks <- NULL
  last_fit <- NULL
  last_k <- NA
  
  for (i in L:(N - M)) {
    past <- 1:i
    oos <- (i + 1):(i + M)
    
    if (is_null(last_fit) || last_k > k_threshold) {
      i_refit <- i
      refits <- c(refits, i)
      # last_fit <-
      loglik <- last_fit$draws("trial_log_lik", format = "matrix")
    }
  }
}


lfo_pointwise_log_lik <- function(res) {
  # more stable than log(sum(exp(x))) 
  log_sum_exp <- function(x) {
    max_x <- max(x)  
    max_x + log(sum(exp(x - max_x)))
  }
  
  # more stable than log(mean(exp(x)))
  log_mean_exp <- function(x) {
    log_sum_exp(x) - log(length(x))
  } 
  
  res |> 
    spread_rvars(oos_log_lik[n, i]) |>
    filter(min(oos_log_lik) < 0) |> 
    group_by(n, i) |> 
    mutate(E_log_lik = log_mean_exp(oos_log_lik))  
}

lfo_log_lik <- function(res) {
  res |> 
    lfo_pointwise_log_lik() |> 
    group_by(n) |> 
    summarize(
      E_log_lik = sum(E_log_lik), 
      log_ratio = posterior::rvar_sum(oos_log_lik), 
      .groups = "drop"
    ) |> 
    mutate(
      psis_obj = map(log_ratio, \(lr) suppressWarnings(loo::psis(posterior::draws_of(lr)))),
      k = map_dbl(psis_obj, loo::pareto_k_values)
    )  
}
