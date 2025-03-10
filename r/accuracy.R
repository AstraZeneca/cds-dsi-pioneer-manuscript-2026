get_trial_loo <- function(res, log_lik_var = "trial_log_lik", moment_match = TRUE, ...) {
  res$loo(log_lik_var, moment_match = moment_match, save_psis = TRUE, ...)
}

add_stacked_results <- function(res_data, stacking_weights, ..., by = c("model_type", "trial")) {
  inner_join(res_data, stacking_weights, by = by) |>  
    mutate(weight = as.numeric(weight)) %>%
    bind_rows(
      filter(., !is.na(weight)) |> 
        group_by(fit_type, trial) |> 
        summarize(model_type = "stacked", across(c(...), \(res) stack_draws(res, weight))) 
    )
}

get_trial_c_index <- function(res, analysis_data = NULL) {
  if (!is_null(analysis_data)) { 
    res <- recover_types(res, select(analysis_data, trial))
  }
  
  spread_rvars(res, trial_c_index[trial])
}

get_patient_pointwise_loo <- function(model_loo, stan_data) {
  as_tibble(stan_data[c("patient", "patient_trial")]) |> 
    rename(trial = patient_trial) |> 
    mutate(
      pareto_k_influence = loo::pareto_k_influence_values(model_loo), imputed = row_number() %in% stan_data$imputed_patients,
      elpd_loo = loo::pointwise(model_loo, "elpd_loo")
    )
}

simplex_allocate <- function(simplex, total) {
  # Input validation
  if (abs(sum(simplex) - 1) > 1e-5) {
    stop("Input vector must sum to 1")
  }
  if (total %% 1 != 0) {
    stop("total must be an integer")
  }
  
  # Initial allocation using floor after multiplication
  raw_allocation <- simplex * total
  initial_allocation <- floor(raw_allocation)
  
  # Calculate remaining amount to distribute
  remainder <- total - sum(initial_allocation)
  
  if (remainder > 0) {
    # Get fractional parts
    fractional_parts <- raw_allocation - initial_allocation
    # Get indices that would sort in descending order
    sorted_indices <- order(fractional_parts, decreasing = TRUE)
    
    # Only distribute up to the remainder amount
    result <- initial_allocation
    if (remainder > 0) {
      result[sorted_indices[1:remainder]] <- result[sorted_indices[1:remainder]] + 1
    }
    return(result)
  } else {
    return(initial_allocation)
  }
}

stack_draws <- function(rvs, simplex) {
  map2(rvs, simplex_allocate(simplex, ndraws(rvs[1])), \(rv, n) resample_draws(rv, ndraws = n)) |> 
    map(\(d) as.vector(draws_of(d))) |> 
    purrr::flatten_dbl() |> 
    rvar()
}

get_loo_admin_brier_score <- function(res, loo_obj) {
  res |> 
    spread_draws(trial_admin_brier_score[i, t]) |> 
    ungroup() |> 
    select(.draw, i, t, trial_admin_brier_score) |> 
    pivot_wider(id_cols = c(.draw, t), names_from = i, values_from = trial_admin_brier_score) |> 
    select(!.draw) |> 
    nest(draws_matrix = !t) |> 
    transmute(t, mean_brier_score = map_dbl(draws_matrix, \(m) sum(loo::E_loo(as.matrix(m), loo_obj$psis_object, type = "mean")$value)))
}
