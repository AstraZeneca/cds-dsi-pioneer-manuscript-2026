library(testthat)
library(here)
library(stringr)

# ============================================================================
# Helpers
# ============================================================================

run_multistate_flags <- function(data) {
  test_stan_function(
    here("tests", "testthat", "stan", "test_multistate_flags_all.stan"),
    data
  )
}

# ============================================================================
# B1: compute_ms_time_scale_flags
# ============================================================================

test_that("compute_ms_time_scale_flags: correct routing for all flag combos", {
  enable_vals    <- c(0L, 1L)
  time_scale_vals <- c(-1L, 0L, 1L, 2L, 3L)

  combos <- expand.grid(
    enable_ms_12    = enable_vals,
    ms_time_scale_12 = time_scale_vals
  )

  n_combos <- nrow(combos)

  data <- list(
    n_combos_b1             = n_combos,
    enable_ms_12_combos     = combos$enable_ms_12,
    ms_time_scale_12_combos = combos$ms_time_scale_12,
    # B2 minimal (not tested here — set n_combos_b2=0 via n_levels_b2 hack)
    # We still need valid B2/B3/B4 inputs; use n_combos=1 trivial case
    n_levels_b2                        = 1L,
    n_combos_b2                        = 1L,
    n_forecast_groups_b2               = array(2L, dim = 1L),
    enable_ms_level_baseline_combos    = array(1L, dim = c(1L, 1L)),
    n_enabled_b3                       = 2L,
    n_gp_b3                            = 0L,
    n_combos_b3                        = 1L,
    enable_ms_01_combos_b3             = array(1L, dim = 1L),
    enable_ms_02_combos_b3             = array(0L, dim = 1L),
    need_12_s_gp_combos_b3            = array(0L, dim = 1L),
    need_12_t_gp_combos_b3            = array(0L, dim = 1L),
    enable_ms_03_combos_b3             = array(0L, dim = 1L),
    enable_ms_32_combos_b3             = array(0L, dim = 1L),
    n_combos_b4                        = 1L,
    n_visits_b4                        = 10L,
    enable_ms_visit_gated_01_combos_b4 = array(0L, dim = 1L),
    enable_ms_visit_gated_latent_01_combos_b4 = array(0L, dim = 1L)
  )

  fit <- run_multistate_flags(data)
  d   <- posterior::as_draws_df(fit$draws())

  for (ci in seq_len(n_combos)) {
    en  <- combos$enable_ms_12[ci]
    sc  <- combos$ms_time_scale_12[ci]
    lbl <- sprintf("enable=%d scale=%d", en, sc)

    if (sc < 0 || sc > 2) {
      # Invalid: is_valid == 0
      expect_equal(get_stan_val(d, "out_b1_is_valid", ci), 0L,
        label = paste(lbl, "is_valid"))
    } else {
      expect_equal(get_stan_val(d, "out_b1_is_valid", ci), 1L,
        label = paste(lbl, "is_valid"))

      exp_need_s <- as.integer(en == 1L && sc %in% c(1L, 2L))
      exp_need_t <- as.integer(en == 1L && sc %in% c(0L, 2L))
      exp_t_int  <- as.integer(exp_need_t == 1L && exp_need_s == 0L)

      expect_equal(get_stan_val(d, "out_b1_need_s", ci), exp_need_s,
        label = paste(lbl, "need_s"))
      expect_equal(get_stan_val(d, "out_b1_need_t", ci), exp_need_t,
        label = paste(lbl, "need_t"))
      expect_equal(get_stan_val(d, "out_b1_t_int",  ci), exp_t_int,
        label = paste(lbl, "t_int"))

      # Invariant: when disabled, both GPs off
      if (en == 0L) {
        expect_equal(get_stan_val(d, "out_b1_need_s", ci), 0L,
          label = paste(lbl, "need_s disabled"))
        expect_equal(get_stan_val(d, "out_b1_need_t", ci), 0L,
          label = paste(lbl, "need_t disabled"))
      }

      # Extended: both GPs on, t_int off
      if (en == 1L && sc == 2L) {
        expect_equal(get_stan_val(d, "out_b1_need_s", ci), 1L,
          label = paste(lbl, "extended need_s"))
        expect_equal(get_stan_val(d, "out_b1_need_t", ci), 1L,
          label = paste(lbl, "extended need_t"))
        expect_equal(get_stan_val(d, "out_b1_t_int",  ci), 0L,
          label = paste(lbl, "extended t_int"))
      }
    }
  }
})

# ============================================================================
# B2: compute_ms_level_baseline_flags
# ============================================================================

test_that("compute_ms_level_baseline_flags: correct GP routing for all hazard flag combos", {
  n_levels_b2 <- 2L
  # 2 groups at each level
  n_forecast_groups_b2 <- c(2L, 3L)

  # Valid combos: all permutations of {0,1,2,3} x {0,1,2,3}
  valid_grid <- expand.grid(lv1 = 0:3, lv2 = 0:3)
  n_valid <- nrow(valid_grid)

  # Invalid combos: at least one level has -1 or 4
  invalid_grid <- expand.grid(lv1 = c(-1L, 4L), lv2 = c(0L, 2L))
  n_invalid    <- nrow(invalid_grid)

  all_combos <- rbind(valid_grid, invalid_grid)
  n_combos   <- nrow(all_combos)

  data <- list(
    n_combos_b1             = 1L,
    enable_ms_12_combos     = array(1L, dim = 1L),
    ms_time_scale_12_combos = array(1L, dim = 1L),
    n_levels_b2                        = n_levels_b2,
    n_combos_b2                        = n_combos,
    n_forecast_groups_b2               = n_forecast_groups_b2,
    enable_ms_level_baseline_combos    = as.matrix(all_combos),
    n_enabled_b3                       = 2L,
    n_gp_b3                            = 0L,
    n_combos_b3                        = 1L,
    enable_ms_01_combos_b3             = array(1L, dim = 1L),
    enable_ms_02_combos_b3             = array(0L, dim = 1L),
    need_12_s_gp_combos_b3            = array(0L, dim = 1L),
    need_12_t_gp_combos_b3            = array(0L, dim = 1L),
    enable_ms_03_combos_b3             = array(0L, dim = 1L),
    enable_ms_32_combos_b3             = array(0L, dim = 1L),
    n_combos_b4                        = 1L,
    n_visits_b4                        = 10L,
    enable_ms_visit_gated_01_combos_b4 = array(0L, dim = 1L),
    enable_ms_visit_gated_latent_01_combos_b4 = array(0L, dim = 1L)
  )

  fit <- run_multistate_flags(data)
  d   <- posterior::as_draws_df(fit$draws())

  for (ci in seq_len(n_combos)) {
    flags <- as.integer(all_combos[ci, ])  # length n_levels_b2
    lbl   <- sprintf("combo %d [%s]", ci, paste(flags, collapse = ","))

    if (any(flags < 0 | flags > 3)) {
      expect_equal(get_stan_val(d, "out_b2_is_valid", ci), 0L,
        label = paste(lbl, "is_valid"))
    } else {
      expect_equal(get_stan_val(d, "out_b2_is_valid", ci), 1L,
        label = paste(lbl, "is_valid"))

      # any_re: max >= 2
      exp_any_re <- as.integer(max(flags) >= 2L)
      expect_equal(get_stan_val(d, "out_b2_any_re", ci), exp_any_re,
        label = paste(lbl, "any_re"))

      # is_gp per level
      exp_is_gp <- as.integer(flags == 3L)
      for (lv in seq_len(n_levels_b2)) {
        expect_equal(get_stan_val(d, "out_b2_is_gp", ci, lv), exp_is_gp[lv],
          label = paste(lbl, sprintf("is_gp[%d]", lv)))
      }

      # n_gp_groups via R oracle
      exp_n_gp <- r_compute_n_enabled_groups(n_forecast_groups_b2, exp_is_gp)
      expect_equal(get_stan_val(d, "out_b2_n_gp", ci), exp_n_gp,
        label = paste(lbl, "n_gp"))

      # n_gp_groups == 0 when no level is 3
      if (!any(flags == 3L)) {
        expect_equal(get_stan_val(d, "out_b2_n_gp", ci), 0L,
          label = paste(lbl, "n_gp when no GP levels"))
      }

      # gp_pos via R oracle
      exp_gp_pos <- r_create_enabled_pos(n_forecast_groups_b2, exp_is_gp)
      for (k in seq_len(n_levels_b2 + 1L)) {
        expect_equal(get_stan_val(d, "out_b2_gp_pos", ci, k), exp_gp_pos[k],
          label = paste(lbl, sprintf("gp_pos[%d]", k)))
      }

      # enabled_pos (any flag >= 1)
      exp_is_enabled <- as.integer(flags >= 1L)
      exp_en_pos <- r_create_enabled_pos(n_forecast_groups_b2, exp_is_enabled)
      for (k in seq_len(n_levels_b2 + 1L)) {
        expect_equal(get_stan_val(d, "out_b2_enabled_pos", ci, k), exp_en_pos[k],
          label = paste(lbl, sprintf("enabled_pos[%d]", k)))
      }
    }
  }
})

# ============================================================================
# B3: compute_ms_transition_group_counts
# ============================================================================

test_that("compute_ms_transition_group_counts: correct group counts for all transition combos", {
  n_enabled <- 5L
  n_gp      <- 2L

  # All 64 combos of 6 binary flags
  all_b3 <- expand.grid(
    en01  = 0:1, en02  = 0:1, s_gp  = 0:1,
    t_gp  = 0:1, en03  = 0:1, en32  = 0:1
  )

  # Separate valid (en32=0 | en03=1) from invalid
  valid_b3   <- all_b3[!(all_b3$en32 == 1L & all_b3$en03 == 0L), ]
  invalid_b3 <- all_b3[  all_b3$en32 == 1L & all_b3$en03 == 0L, , drop = FALSE]
  combos_b3  <- rbind(valid_b3, invalid_b3)
  n_combos   <- nrow(combos_b3)

  data <- list(
    n_combos_b1             = 1L,
    enable_ms_12_combos     = array(1L, dim = 1L),
    ms_time_scale_12_combos = array(1L, dim = 1L),
    n_levels_b2                        = 1L,
    n_combos_b2                        = 1L,
    n_forecast_groups_b2               = array(2L, dim = 1L),
    enable_ms_level_baseline_combos    = array(1L, dim = c(1L, 1L)),
    n_enabled_b3                       = n_enabled,
    n_gp_b3                            = n_gp,
    n_combos_b3                        = n_combos,
    enable_ms_01_combos_b3             = combos_b3$en01,
    enable_ms_02_combos_b3             = combos_b3$en02,
    need_12_s_gp_combos_b3            = combos_b3$s_gp,
    need_12_t_gp_combos_b3            = combos_b3$t_gp,
    enable_ms_03_combos_b3             = combos_b3$en03,
    enable_ms_32_combos_b3             = combos_b3$en32,
    n_combos_b4                        = 1L,
    n_visits_b4                        = 10L,
    enable_ms_visit_gated_01_combos_b4 = array(0L, dim = 1L),
    enable_ms_visit_gated_latent_01_combos_b4 = array(0L, dim = 1L)
  )

  fit <- run_multistate_flags(data)
  d   <- posterior::as_draws_df(fit$draws())

  n_valid_combos   <- nrow(valid_b3)
  n_invalid_combos <- nrow(invalid_b3)

  # Check valid combos
  for (ci in seq_len(n_valid_combos)) {
    row <- valid_b3[ci, ]
    lbl <- sprintf("combo %d [en01=%d en02=%d sgp=%d tgp=%d en03=%d en32=%d]",
      ci, row$en01, row$en02, row$s_gp, row$t_gp, row$en03, row$en32)

    expect_equal(get_stan_val(d, "out_b3_is_valid", ci), 1L,
      label = paste(lbl, "is_valid"))

    transitions <- list(
      list(flag = row$en01,  en_name = "out_b3_n_en_01",  gp_name = "out_b3_n_gp_01"),
      list(flag = row$en02,  en_name = "out_b3_n_en_02",  gp_name = "out_b3_n_gp_02"),
      list(flag = row$s_gp,  en_name = "out_b3_n_en_12s", gp_name = "out_b3_n_gp_12s"),
      list(flag = row$t_gp,  en_name = "out_b3_n_en_12t", gp_name = "out_b3_n_gp_12t"),
      list(flag = row$en03,  en_name = "out_b3_n_en_03",  gp_name = "out_b3_n_gp_03"),
      list(flag = row$en32,  en_name = "out_b3_n_en_32",  gp_name = "out_b3_n_gp_32")
    )

    for (tr in transitions) {
      exp_en <- if (tr$flag == 1L) n_enabled else 0L
      exp_gp <- if (tr$flag == 1L) n_gp      else 0L

      expect_equal(get_stan_val(d, tr$en_name, ci), exp_en,
        label = paste(lbl, tr$en_name))
      expect_equal(get_stan_val(d, tr$gp_name, ci), exp_gp,
        label = paste(lbl, tr$gp_name))

      # Invariant: n_gp <= n_enabled
      expect_true(exp_gp <= exp_en,
        label = paste(lbl, "n_gp <= n_enabled for", tr$en_name))
    }
  }

  # Check invalid combos (en32=1, en03=0) — placed after valid ones
  for (j in seq_len(n_invalid_combos)) {
    ci <- n_valid_combos + j
    expect_equal(get_stan_val(d, "out_b3_is_valid", ci), 0L,
      label = sprintf("invalid combo %d is_valid", ci))
  }
})

# ============================================================================
# B4: compute_ms_obs_psa_covar_size
# ============================================================================

test_that("compute_ms_obs_psa_covar_size: size is n_visits only for visit_gated=1, latent=0", {
  n_visits <- 42L

  combos_b4 <- expand.grid(
    vg  = 0:1,
    lat = 0:1
  )
  n_combos <- nrow(combos_b4)

  data <- list(
    n_combos_b1             = 1L,
    enable_ms_12_combos     = array(1L, dim = 1L),
    ms_time_scale_12_combos = array(1L, dim = 1L),
    n_levels_b2                        = 1L,
    n_combos_b2                        = 1L,
    n_forecast_groups_b2               = array(2L, dim = 1L),
    enable_ms_level_baseline_combos    = array(1L, dim = c(1L, 1L)),
    n_enabled_b3                       = 2L,
    n_gp_b3                            = 0L,
    n_combos_b3                        = 1L,
    enable_ms_01_combos_b3             = array(1L, dim = 1L),
    enable_ms_02_combos_b3             = array(0L, dim = 1L),
    need_12_s_gp_combos_b3            = array(0L, dim = 1L),
    need_12_t_gp_combos_b3            = array(0L, dim = 1L),
    enable_ms_03_combos_b3             = array(0L, dim = 1L),
    enable_ms_32_combos_b3             = array(0L, dim = 1L),
    n_combos_b4                        = n_combos,
    n_visits_b4                        = n_visits,
    enable_ms_visit_gated_01_combos_b4 = combos_b4$vg,
    enable_ms_visit_gated_latent_01_combos_b4 = combos_b4$lat
  )

  fit <- run_multistate_flags(data)
  d   <- posterior::as_draws_df(fit$draws())

  for (ci in seq_len(n_combos)) {
    vg  <- combos_b4$vg[ci]
    lat <- combos_b4$lat[ci]
    lbl <- sprintf("vg=%d lat=%d", vg, lat)

    exp_size <- if (vg == 1L && lat == 0L) n_visits else 0L
    expect_equal(get_stan_val(d, "out_b4_size", ci), exp_size,
      label = paste(lbl, "size"))
  }
})
