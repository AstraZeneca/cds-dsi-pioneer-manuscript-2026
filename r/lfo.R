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