library(testthat)
library(here)
library(posterior)
library(tidyverse)

r_full_model_grid_flags <- function(pop_pn, pat_pn, states_grid,
                                     ms_tv_cov, n_tv, ms_01, visit_gated,
                                     visit_gated_latent, tv_02_cov, tv_03_cov,
                                     has_inline_tv_covar) {
  any_pn <- as.integer(pop_pn || pat_pn)
  needs_inline <- as.integer(
    has_inline_tv_covar &&
    !any_pn &&
    ms_tv_cov && n_tv > 0 &&
    (tv_02_cov || tv_03_cov || (visit_gated && visit_gated_latent))
  )
  need_grid <- as.integer(
    any_pn || states_grid ||
    (ms_tv_cov && n_tv > 0 && !needs_inline &&
     ((ms_01 && !visit_gated) ||
      (visit_gated && visit_gated_latent) ||
      tv_02_cov || tv_03_cov))
  )
  list(any_pn = any_pn, need_grid = need_grid, needs_inline = needs_inline)
}

combos <- expand.grid(
  pop_pn              = 0:1,
  pat_pn              = 0:1,
  states_grid         = 0:1,
  ms_tv_cov           = 0:1,
  n_tv                = c(0L, 1L, 3L),
  ms_01               = 0:1,
  visit_gated         = 0:1,
  visit_gated_latent  = 0:1,
  tv_02_cov           = 0:1,
  tv_03_cov           = 0:1,
  has_inline_tv_covar = 0:1
) |> as.data.frame()

fit <- test_stan_function(
  here("tests", "testthat", "stan", "test_full_model_flags_all.stan"),
  list(
    n_combos                    = nrow(combos),
    enable_pop_pn               = combos$pop_pn,
    enable_patient_pn           = combos$pat_pn,
    enable_states_grid          = combos$states_grid,
    enable_ms_tv_cov            = combos$ms_tv_cov,
    n_tv_covar                  = combos$n_tv,
    enable_ms_01                = combos$ms_01,
    enable_visit_gated          = combos$visit_gated,
    enable_visit_gated_latent   = combos$visit_gated_latent,
    enable_02_tv_cov            = combos$tv_02_cov,
    enable_03_tv_cov            = combos$tv_03_cov,
    has_inline_tv_covar         = combos$has_inline_tv_covar
  )
)
d <- posterior::as_draws_df(fit$draws())

test_that("compute_full_model_grid_flags: all combos correct", {
  for (c in seq_len(nrow(combos))) {
    r   <- combos[c, ]
    exp <- r_full_model_grid_flags(r$pop_pn, r$pat_pn, r$states_grid,
                                    r$ms_tv_cov, r$n_tv, r$ms_01,
                                    r$visit_gated, r$visit_gated_latent,
                                    r$tv_02_cov, r$tv_03_cov,
                                    r$has_inline_tv_covar)
    expect_equal(get_stan_val(d, "out_any_process_noise", c),     exp$any_pn,
                 label = sprintf("any_pn[%d]", c))
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), exp$need_grid,
                 label = sprintf("need_grid[%d]", c))
    expect_equal(get_stan_val(d, "out_ms_needs_inline_burden", c),   exp$needs_inline,
                 label = sprintf("needs_inline[%d]", c))
  }
})

test_that("compute_full_model_grid_flags: need_grid=1 when any process noise", {
  pn_combos <- which(combos$pop_pn == 1L | combos$pat_pn == 1L)
  for (c in pn_combos)
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), 1L,
                 label = sprintf("need_grid pn[%d]", c))
})

test_that("compute_full_model_grid_flags: needs_inline=0 when process noise on", {
  pn_combos <- which(combos$pop_pn == 1L | combos$pat_pn == 1L)
  for (c in pn_combos)
    expect_equal(get_stan_val(d, "out_ms_needs_inline_burden", c), 0L,
                 label = sprintf("needs_inline pn[%d]", c))
})

test_that("compute_full_model_grid_flags: needs_inline=0 when caller has no inline path", {
  # has_inline_tv_covar = 0 must always disable the inline path, regardless of
  # what other flags say. Tumor models depend on this — they have no
  # inline TV-covariate populator, so going down the inline path leaves
  # ms_time_varying_covar_01 unpopulated (with downstream out-of-range crashes).
  no_inline_combos <- which(combos$has_inline_tv_covar == 0L)
  for (c in no_inline_combos)
    expect_equal(get_stan_val(d, "out_ms_needs_inline_burden", c), 0L,
                 label = sprintf("needs_inline (no inline path)[%d]", c))
})

test_that("compute_full_model_grid_flags: latent visit-gated forces grid when no inline path", {
  # The case that was broken before adding has_inline_tv_covar: tumor model
  # with latent visit-gating + TV cov on. Without the gate, needs_inline
  # incorrectly fired and need_grid came back 0 → empty states_full_grid.
  vg_latent_combos <- which(
    combos$has_inline_tv_covar == 0L &
    combos$visit_gated == 1L &
    combos$visit_gated_latent == 1L &
    combos$ms_tv_cov == 1L &
    combos$n_tv > 0L &
    combos$pop_pn == 0L & combos$pat_pn == 0L
  )
  for (c in vg_latent_combos) {
    expect_equal(get_stan_val(d, "out_ms_needs_inline_burden", c), 0L,
                 label = sprintf("needs_inline tumor latent[%d]", c))
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), 1L,
                 label = sprintf("need_grid tumor latent[%d]", c))
  }
})

test_that("compute_full_model_grid_flags: need_grid=0 when all gates off", {
  off_combos <- which(vapply(seq_len(nrow(combos)), function(c) {
    r <- combos[c, ]
    exp <- r_full_model_grid_flags(r$pop_pn, r$pat_pn, r$states_grid,
                                   r$ms_tv_cov, r$n_tv, r$ms_01,
                                   r$visit_gated, r$visit_gated_latent,
                                   r$tv_02_cov, r$tv_03_cov,
                                   r$has_inline_tv_covar)
    exp$need_grid == 0L
  }, logical(1)))
  for (c in off_combos)
    expect_equal(get_stan_val(d, "out_need_states_full_grid", c), 0L,
                 label = sprintf("need_grid off[%d]", c))
})
