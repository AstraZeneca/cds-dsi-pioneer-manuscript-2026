# This file runs all test cases in a single Stan model call, using a rectangular output array.
#
library(here)
library(testthat)

# Define RECIST constants for readability
CR <- 1 # Complete Response
PR <- 2 # Partial Response
SD <- 3 # Stable Disease
PD <- 4 # Progressive Disease

# List of test cases: each is a list with sld_trajectory, n_screening, pre_nadir, expected, expected_with_nadir, description
recist_cases <- list(
  # Edge: SLD exactly at 20% increase and exactly 5mm above nadir (PD threshold)
  list(
    sld_trajectory = c(100, 80, 96),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, PD),
    expected_with_nadir = c(PR, PD),
    desc = "PD: exactly 20% and 5mm above nadir"
  ),
  # Edge: SLD exactly at 30% decrease from baseline (PR threshold)
  list(
    sld_trajectory = c(100, 70),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR),
    expected_with_nadir = c(PR),
    desc = "PR: exactly 30% decrease from baseline"
  ),
  # Edge: SLD just below 20% increase or just below 5mm above nadir (should not be PD)
  list(
    sld_trajectory = c(100, 80, 95.9),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "Not PD: just below 20% or 5mm"
  ),
  # Edge: SLD just above 20% increase but <5mm above nadir (should not be PD)
  list(
    sld_trajectory = c(100, 80, 84.1),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "Not PD: >20% but <5mm above nadir"
  ),
  # Edge: SLD just above 5mm but <20% increase from nadir (should not be PD)
  list(
    sld_trajectory = c(100, 80, 85.1),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "Not PD: >5mm but <20% above nadir"
  ),
  # Edge: SLD returns to baseline after PR or CR (should revert to SD, not PR/CR)
  list(
    sld_trajectory = c(100, 60, 100),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "Return to baseline after PR"
  ),
  # Edge: SLD = 0 after PD (should be CR)
  list(
    sld_trajectory = c(100, 130, 0),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PD, CR),
    expected_with_nadir = c(PD, CR),
    desc = "CR after PD"
  ),
  # Edge: SLD = baseline at all timepoints (should always be SD)
  list(
    sld_trajectory = c(100, 100, 100, 100),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD, SD, SD),
    expected_with_nadir = c(SD, SD, SD),
    desc = "Always SD: SLD = baseline"
  ),
  # Edge: SLD fluctuates around nadir and baseline, never crossing thresholds
  list(
    sld_trajectory = c(100, 95, 105, 98, 102),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD, SD, SD, SD),
    expected_with_nadir = c(SD, SD, SD, SD),
    desc = "Fluctuating SLD, always SD"
  ),
  # Edge: SLD = 0 at baseline, then increases (should be CR, then SD/PD as appropriate)
  list(
    sld_trajectory = c(0, 0, 10, 20),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(CR, SD, SD),
    expected_with_nadir = c(CR, SD, SD),
    desc = "Baseline CR, then increase"
  ),
  # Basic RECIST classifications
  list(
    sld_trajectory = c(100, 0),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(CR),
    expected_with_nadir = c(CR),
    desc = "Complete response"
  ),
  list(
    sld_trajectory = c(100, 70),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR),
    expected_with_nadir = c(PR),
    desc = "Partial response: exactly 30%"
  ),
  list(
    sld_trajectory = c(100, 60),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR),
    expected_with_nadir = c(PR),
    desc = "Partial response: >30%"
  ),
  list(
    sld_trajectory = c(100, 80),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD),
    expected_with_nadir = c(SD),
    desc = "Stable disease: 20% decrease"
  ),
  list(
    sld_trajectory = c(100, 100),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD),
    expected_with_nadir = c(SD),
    desc = "Stable disease: no change"
  ),
  list(
    sld_trajectory = c(100, 50, 65),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, PD),
    expected_with_nadir = c(PR, PD),
    desc = "Progressive disease"
  ),

  # Boundary conditions
  list(
    sld_trajectory = c(100, 70.1),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD),
    expected_with_nadir = c(SD),
    desc = "Just below PR threshold"
  ),
  list(
    sld_trajectory = c(100, 69.9),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR),
    expected_with_nadir = c(PR),
    desc = "Just above PR threshold"
  ),
  list(
    sld_trajectory = c(20, 8, 10),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "25% increase but <5mm"
  ),
  list(
    sld_trajectory = c(100, 40, 48),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, PD),
    expected_with_nadir = c(PR, PD),
    desc = "20% increase and >5mm"
  ),

  # Nadir tracking
  list(
    sld_trajectory = c(100, 80, 60, 70, 85),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD, PR, SD, PD),
    expected_with_nadir = c(SD, PR, SD, PD),
    desc = "Nadir tracking"
  ),
  list(
    sld_trajectory = c(100, 70, 90, 60, 80),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, SD, PR, SD),
    expected_with_nadir = c(PR, SD, PR, SD),
    desc = "Multiple decreases"
  ),

  # Multiple screening visits
  list(
    sld_trajectory = c(120, 100, 70),
    n_screening = 2,
    pre_nadir = 0,
    expected = c(PR),
    expected_with_nadir = c(PR),
    desc = "2 screening visits"
  ),
  list(
    sld_trajectory = c(110, 105, 100, 60, 80),
    n_screening = 3,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "3 screening visits"
  ),
  list(
    sld_trajectory = c(100, 95, 70, 85),
    n_screening = 2,
    pre_nadir = 0,
    expected = c(PR, SD),
    expected_with_nadir = c(PR, SD),
    desc = "Output excludes screening"
  ),

  # Clinical scenarios
  list(
    sld_trajectory = c(100, 65, 70, 68, 72),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, PR, PR, PR),
    expected_with_nadir = c(PR, PR, PR, PR),
    desc = "Sustained response"
  ),
  list(
    sld_trajectory = c(100, 60, 65, 85),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(PR, PR, PD),
    expected_with_nadir = c(PR, PR, PD),
    desc = "Response then progression"
  ),
  list(
    sld_trajectory = c(100, 105, 110, 125),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD, SD, PD),
    expected_with_nadir = c(SD, SD, PD),
    desc = "Gradual progression"
  ),
  list(
    sld_trajectory = c(80, 0, 0, 25),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(CR, CR, SD),
    expected_with_nadir = c(CR, CR, SD),
    desc = "CR with recurrence"
  ),

  # Pre-nadir function signature
  list(
    sld_trajectory = c(100, 80, 60, 85),
    n_screening = 1,
    pre_nadir = 70,
    expected = c(SD, PR, SD),
    expected_with_nadir = c(SD, PR, SD),
    desc = "With pre-nadir=70"
  ),
  list(
    sld_trajectory = c(100, 80, 60, 85),
    n_screening = 1,
    pre_nadir = 0,
    expected = c(SD, PR, SD),
    expected_with_nadir = c(SD, PR, SD),
    desc = "Auto nadir"
  )
)

max_len <- max(sapply(recist_cases, function(x) length(x$sld_trajectory)))
n_cases <- length(recist_cases)

stan_data <- list(
  n_cases = n_cases,
  max_len = max_len,
  n_trajectory = sapply(recist_cases, function(x) length(x$sld_trajectory)),
  n_screening = sapply(recist_cases, function(x) x$n_screening),
  sld_trajectory = t(sapply(recist_cases, function(x) {
    v <- rep(0, max_len)
    v[1:length(x$sld_trajectory)] <- x$sld_trajectory
    v
  })),
  pre_nadir = sapply(recist_cases, function(x) x$pre_nadir)
)

fit <- test_stan_function(
  here::here(
    "tests",
    "testthat",
    "stan",
    "test_calculate_target_recist_all.stan"
  ),
  stan_data
)


# Stan outputs as flat matrix: [draw, variable], variable names are result[case,time]
result_matrix <- fit$draws("result", format = "matrix")
result_with_nadir_matrix <- fit$draws("result_with_nadir", format = "matrix")
# Use first draw, reshape to [n_cases, max_len]
n_cases <- stan_data$n_cases
max_len <- stan_data$max_len
result <- matrix(
  result_matrix[1, ],
  nrow = n_cases,
  ncol = max_len,
  byrow = TRUE
)
result_with_nadir <- matrix(
  result_with_nadir_matrix[1, ],
  nrow = n_cases,
  ncol = max_len,
  byrow = TRUE
)

for (i in seq_along(recist_cases)) {
  n_treat <- stan_data$n_trajectory[i] - stan_data$n_screening[i]
  # Only check up to the available columns in the output
  check_len <- min(n_treat, ncol(result))
  diagnosed <- FALSE
  test_that(sprintf("calculate_target_recist: %s", recist_cases[[i]]$desc), {
    r_out <- as.integer(result[i, 1:check_len])
    stan_out <- as.integer(result_with_nadir[i, 1:check_len])
    if (!all(r_out == stan_out) && !diagnosed) {
      diagnosed <<- TRUE
      cat("\n--- DIAGNOSTIC OUTPUT FOR FIRST FAILING CASE ---\n")
      cat("Case:", i, "\n")
      cat(
        "SLD:",
        str_c(recist_cases[[i]]$sld_trajectory, collapse = ", "),
        "\n"
      )
      cat("R output:", str_c(r_out, collapse = ", "), "\n")
      cat("Stan output:", str_c(stan_out, collapse = ", "), "\n")
      # Step-by-step for R
      sld <- recist_cases[[i]]$sld_trajectory
      nadir <- cummin(sld)
      baseline <- sld[1]
      cat("Step-by-step (R):\n")
      for (j in 2:length(sld)) {
        current <- sld[j]
        curr_nadir <- nadir[j]
        cr <- current == 0
        pd <- curr_nadir > 0 &&
          current > curr_nadir &&
          (current - curr_nadir) / curr_nadir >= 0.2 &&
          (current - curr_nadir) >= 5
        pr <- (current - baseline) / baseline <= -0.3
        cat(sprintf(
          "  t=%d: SLD=%.2f, nadir=%.2f, CR=%s, PD=%s, PR=%s\n",
          j,
          current,
          curr_nadir,
          cr,
          pd,
          pr
        ))
      }
      # Stan logic (should match R, but print for clarity)
      cat("Step-by-step (Stan logic):\n")
      nadir_stan <- sld[1]
      for (j in 2:length(sld)) {
        nadir_stan <- min(nadir_stan, sld[j])
        current <- sld[j]
        cr <- current == 0
        pd <- nadir_stan > 0 &&
          current > nadir_stan &&
          (current - nadir_stan) / nadir_stan >= 0.2 &&
          (current - nadir_stan) >= 5
        pr <- (current - baseline) / baseline <= -0.3
        cat(sprintf(
          "  t=%d: SLD=%.2f, nadir=%.2f, CR=%s, PD=%s, PR=%s\n",
          j,
          current,
          nadir_stan,
          cr,
          pd,
          pr
        ))
      }
      cat("--- END DIAGNOSTIC ---\n\n")
    }
    expect_equal(stan_out, r_out)
  })
}
