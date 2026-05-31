if (!exists("init_project")) source(here::here(".Rprofile"))
init_project()

# The publication pipeline reuses the sclc model code (Stan model,
# multistate logic, initializers, priors, prepare_tumor_stan_data). init_project()
# only sources these when TAR_PROJECT == "sclc", so source them explicitly
# here for TAR_PROJECT == "publication".
source(here::here("r", "sclc", "priors.R"))
source(here::here("r", "multistate.R"))
source(here::here("r", "sclc", "prepare_analysis_data.R"))
source(here::here("r", "sclc", "accuracy.R"))
source(here::here("r", "sclc", "initializers.R"))
source(here::here("r", "publication", "prepare_analysis_data.R"))

publication_data_path <- "/mnt/data/PUBLICATION"
publication_output_path <- tar_path_store() |> fs::path_dir()
publication_artifacts_path <- str_replace(
  publication_output_path,
  "^/mnt/data/analysis-results",
  "/mnt/artifacts/"
)

fs::dir_create(file.path(publication_output_path, "fit"))
fs::dir_create(file.path(publication_artifacts_path, "models"))
fs::dir_create(file.path(publication_artifacts_path, "crew_logs"))

controller_default <- crew_controller_local(name = "default", workers = 4)
controller_fit <- crew_controller_local(name = "fit", workers = 4)
controller_many_samples <- crew_controller_local(name = "many samples", workers = 4)

tar_option_set(
  packages = c(
    "magrittr",
    "tidyverse",
    "rlang",
    "here",
    "targets",
    "cmdstanr",
    "tidybayes",
    "posterior",
    "loo",
    "recipes",
    "arrow"
  ),
  controller = crew_controller_group(
    controller_default,
    controller_fit,
    controller_many_samples
  ),
  resources = tar_resources(crew = tar_resources_crew(controller = "default")),
  memory = "transient",
  garbage_collection = TRUE,
  format = rvar_safe_qs2_format,
  error = "continue",
  seed = 26091468
)

cat("TAR_PROJECT =", Sys.getenv("TAR_PROJECT"), "\n")
cat("Current git branch =", gert::git_branch(), "\n")
cat("User =", Sys.getenv("DOMINO_STARTING_USERNAME"), "\n")
cat("Targets store =", tar_path_store(), "\n")
cat("publication_output_path =", publication_output_path, "\n")

km_quant <- seq(0.2, 0.8, by = 0.05)

pfs_timepoints_pub <- enframe(c(6, 9, 12, 15, 18), name = "n", value = "timepoint")

# SCC covariates: age, sex, ECOG, hgb, LDH, albumin
covar_formula_scc <- ~ age + male + ecog + hgb + ldh_log + albumin
# CRC covariates: ECOG is constant (all zeros) so excluded; otherwise same set
covar_formula_crc <- ~ age + male + hgb + ldh_log + albumin

disease_map <- tibble::tribble(
  ~disease, ~disease_data_path,                    ~covar_formula,    ~trial_re,
  "scc",    publication_data_path,                  covar_formula_scc, FALSE,
  "crc",    file.path(publication_data_path, "crc"), covar_formula_crc, FALSE
)

publication_targets <- list(

  # Track initializer file so changes invalidate the initializer targets
  tar_target(
    initializers_fixed_file,
    "r/sclc/initializers_fixed.R",
    format = "file"
  ),

  # Model (shared across diseases) -----------------------------------------------

  tar_target(
    tumor_ssls_model_file,
    here("stan", "tumor", "sf-ssm-log-space.stan"),
    format = "file"
  ),
  tar_target(
    tumor_ssls_include_files,
    find_stan_includes(tumor_ssls_model_file),
    format = "file"
  ),
  tar_target(
    tumor_ssls_exe_hash,
    build_model_exe_hash(
      tumor_ssls_model_file,
      tumor_ssls_include_files,
      publication_artifacts_path,
      include_paths = c(here("stan"), here("stan", "tumor"))
    ),
    error = "stop",
    cue = tar_cue("always")
  ),

  # Per-disease pipeline ---------------------------------------------------------

  tar_map(
    disease_map,
    names = "disease",

    # Data -----------------------------------------------------------------------

    tar_target(
      target_patient_data_file,
      file.path(disease_data_path, "target", "cooked_patient_data.csv"),
      format = "file"
    ),
    tar_target(
      target_visit_data_file,
      file.path(disease_data_path, "target", "assessment_visit_data.csv"),
      format = "file"
    ),
    tar_target(
      historical_patient_data_file,
      file.path(disease_data_path, "historical", "cooked_patient_data.csv"),
      format = "file"
    ),
    tar_target(
      historical_visit_data_file,
      file.path(disease_data_path, "historical", "assessment_visit_data.csv"),
      format = "file"
    ),

    tar_target(
      target_patient_data,
      read_csv(target_patient_data_file, col_types = cols(studyid = "c", usubjid = "c"), show_col_types = FALSE) |>
        mutate(across(where(is.character), as_factor))
    ),
    tar_target(
      historical_patient_data,
      read_csv(historical_patient_data_file, col_types = cols(studyid = "c", usubjid = "c"), show_col_types = FALSE) |>
        mutate(across(where(is.character), as_factor))
    ),

    tar_target(
      target_visit_data,
      read_csv(target_visit_data_file, col_types = cols(studyid = "c", usubjid = "c"), show_col_types = FALSE) |>
        filter(!is.na(mmsumdiam)) |>
        determine_visit_data_response()
    ),
    tar_target(
      historical_visit_data,
      read_csv(historical_visit_data_file, col_types = cols(studyid = "c", usubjid = "c"), show_col_types = FALSE) |>
        filter(!is.na(mmsumdiam)) |>
        determine_visit_data_response()
    ),

    tar_target(
      all_analysis_data,
      prepare_publication_analysis_data(
        target_patient_data,
        historical_patient_data,
        target_visit_data,
        historical_visit_data
      )
    ),

    # Observed KM curves ---------------------------------------------------------

    tar_target(
      km_trial_pfs,
      get_km_res(all_analysis_data, pfs, right_censored, probs = km_quant)
    ),
    tar_target(
      km_trial_os,
      all_analysis_data |>
        mutate(
          os_time = if_else(death, death_week, patient_max_t),
          os_censored = !death,
          interval_censored = 0L
        ) |>
        get_km_res(os_time, os_censored, probs = km_quant)
    ),

    # Covariates -----------------------------------------------------------------

    tar_target(
      covar_design_matrix,
      prepare_covar_design_matrix(all_analysis_data, covar_formula)
    ),

    tar_target(
      elicited_priors,
      prepare_publication_elicited_priors(
        covar_design_matrix,
        shrink_mean = 1 / 2,
        shrink_sd = 1 / 2
      )
    ),

    # No conditioning subgroups (no pdl1/histology in publication data)
    tar_target(cond_groups, list()),

    tar_target(extend_max_all_t, 200L),

    # Stan data ------------------------------------------------------------------

    tar_target(
      all_stan_data,
      prepare_tumor_stan_data(
        all_analysis_data,
        covar_design_matrix,
        cond_groups,
        km_quant,
        extend_max_all_t = 200L,
        forecast_observation_interval = 6L
      ),
      error = "stop"
    ),

    tar_target(
      default_stan_data_settings,
      lst(
        fit_tumor_data = TRUE,
        fit_multistate_data = TRUE,
        enable_states_full_grid = FALSE,
        sf_rep_T = 20,
        debug = FALSE,

        enable_ms_01 = TRUE,
        enable_ms_02 = TRUE,
        enable_ms_12 = TRUE,
        enable_ms_03 = TRUE,
        enable_ms_32 = TRUE,

        ms_time_scale_12 = 1L,

        # Legacy baseline-hazard mode (per level, 0-4). Kept here as the
        # human-readable config knob; decomposed into the three per-transition x
        # per-level arrays (enable_ms_level_gp / ms_level_intercept_mode /
        # ms_level_intercept_corr_group) at the base_tumor_ssls_stan_data assembly
        # point via decompose_ms_level_baseline_hazard(). corr_group stays all-zero
        # in Phase 1, so the decomposed config reproduces this legacy flag
        # bit-identically.
        enable_ms_level_baseline_hazard = c(trial = 3L, patient = 0L),

        enable_ms_pop_time_varying_cov = TRUE,
        enable_ms_pop_time_invariant_cov = TRUE,
        enable_ms_level_cov = c(trial = FALSE, patient = FALSE),
        # Latent visit-gated 0->1: hazard contributions only at observed visit
        # weeks, but the time-varying covariates (log SLD, log decrease rate,
        # log growth rate) come from the modeled state-space trajectory rather
        # than raw observations. This avoids the log(0) problem that observed
        # mode hits at complete-response visits and aligns the survival
        # likelihood with the actual measurement schedule (lilly_cxcr4 ~6w vs
        # amgen_darbe weekly).
        enable_ms_visit_gated_01 = 1L,
        enable_ms_visit_gated_latent_01 = 1L,
        share_dead_gp_shape = 0L,
        enable_ms_02_time_varying_cov = 0L,
        # 0->3 dropout hazard with patient-level discrimination (added 2026-05-28).
        # The previous fit had no per-patient discrimination on the 0->3 path
        # (only the trial-level GP), so died_off_trial patients were routed
        # to 0->1 too quickly in the spop simulation. Enabling both:
        #   - TI: baseline covariates (age, ECOG, hgb, LDH, albumin, sex)
        #     give static dropout-risk signal.
        #   - TV: latent log SLD / decrease rate / growth rate let dropout
        #     risk track tumor dynamics (e.g. patients on a deteriorating
        #     trajectory may drop out faster).
        # With only 57 dropout events the TV path is identification-limited;
        # priors are kept tight (Normal(0, 0.5)) to avoid overfit. 3->2 TI
        # left off — only 57 events with another competing hazard to model.
        enable_ms_03_time_invariant_cov = 1L,
        enable_ms_03_time_varying_cov = 1L,
        enable_ms_32_time_invariant_cov = 0L,
        enable_ms_12_entry_covar = 0L,
        enable_ms_32_entry_covar = 0L,
        entry_covar_12 = numeric(0),
        entry_covar_32 = numeric(0),

        enable_level_intercept_tr = c(
          trial   = level_intercept_mode[if (trial_re) "re" else "none"],
          patient = level_intercept_mode["re"]
        ),
        enable_level_cov_tr = c(trial = FALSE, patient = FALSE),
        enable_pop_cov_tr = FALSE,
        enable_pop_process_noise_tr = FALSE,
        enable_patient_process_noise_tr = FALSE,
        enable_patient_process_noise_sd_tr = FALSE,
        enable_patient_process_noise_phi_tr = FALSE,

        enable_level_intercept_frac = c(
          trial   = level_intercept_mode[if (trial_re) "re" else "none"],
          patient = level_intercept_mode["re"]
        ),
        enable_level_cov_frac = c(trial = FALSE, patient = FALSE),
        enable_pop_cov_frac = TRUE,

        enable_level_intercept_init = c(
          trial   = level_intercept_mode[if (trial_re) "re" else "none"],
          patient = level_intercept_mode["re"]
        ),
        enable_level_cov_init = c(trial = FALSE, patient = FALSE),
        enable_pop_cov_init = TRUE,

        pfs_timepoints = pfs_timepoints_pub$timepoint,
        n_pfs_timepoints = nrow(pfs_timepoints_pub),

        n_shards = 1L
      )
    ),

    # tumor_priors must be computed *after* default_stan_data_settings so
    # get_tumor_priors() sees enable_ms_visit_gated_01 / enable_ms_02_time_varying_cov
    # (those flags determine the dimension of the time_varying_coef_* hyperprior arrays).
    tar_target(
      tumor_priors,
      get_tumor_priors(
        c(all_stan_data, default_stan_data_settings),
        elicited_priors,
        covar_design_matrix
      )
    ),

    tar_target(
      base_tumor_ssls_stan_data,
      {
        assembled <- all_stan_data |>
          add_tumor_priors(tumor_priors) |>
          c(default_stan_data_settings) |>
          c(derive_ms_fields(all_analysis_data, "full"))
        # Translate the legacy single baseline-hazard flag into the three
        # decomposed per-transition x per-level arrays the Stan model consumes,
        # then drop the legacy key (no longer declared in flags.stan).
        decomposed <- decompose_ms_level_baseline_hazard(
          assembled$enable_ms_level_baseline_hazard
        )
        assembled$enable_ms_level_baseline_hazard <- NULL

        # --- Correlated patient-level frailty on 0->1 and 0->3 (Phase 2) ---
        # Add an RE-NCP patient-level intercept to slots 01 and 03 and place both
        # in correlation group 1 at the patient level. The abundant 0->1
        # progression/censoring history (497 patients) feeds sigma_01; the
        # negatively-learned correlation transmits that evidence to the data-poor
        # 0->3 dropout hazard (57 events), so died_off_trial patients route to
        # dropout (long PFS) instead of a fast 0->1 progression. Slot order:
        # 1=01, 2=02, 3=03, 4=12_s, 5=12_t, 6=32; patient level = last column.
        n_levels_ms <- ncol(decomposed$ms_level_intercept_mode)
        patient_lv <- n_levels_ms
        MS_SLOT_01 <- 1L; MS_SLOT_03 <- 3L
        decomposed$ms_level_intercept_mode[MS_SLOT_01, patient_lv] <- 2L  # RE-NCP
        decomposed$ms_level_intercept_mode[MS_SLOT_03, patient_lv] <- 2L  # RE-NCP
        decomposed$enable_ms_level_gp[MS_SLOT_01, patient_lv] <- 0L
        decomposed$enable_ms_level_gp[MS_SLOT_03, patient_lv] <- 0L
        decomposed$ms_level_intercept_corr_group[MS_SLOT_01, patient_lv] <- 1L
        decomposed$ms_level_intercept_corr_group[MS_SLOT_03, patient_lv] <- 1L

        c(assembled, decomposed)
      }
    ),

    # Fits -----------------------------------------------------------------------

    tar_map(
      tibble(
        type = c("prior", "posterior"),
        fit_data = c(FALSE, TRUE),
        base_name = c("prior_tumor_ssls", "tumor_ssls"),
        iter_sampling = 500,
        iter_warmup = c(300L, 500L),
        chains = 4L
      ),
      names = "type",

      tar_target(
        tumor_ssls_stan_data,
        base_tumor_ssls_stan_data |>
          list_assign(
            fit_tumor_data = fit_data,
            fit_multistate_data = fit_data,
            forecast = TRUE
          )
      ),

      tar_target(
        tumor_ssls_initializer,
        {
          source(initializers_fixed_file)
          create_tumor_ssls_initializer_fixed(tumor_ssls_stan_data)
        }
      ),

      tar_target(
        tumor_ssls_res,
        sample_and_save(
          tumor_ssls_exe_hash$exe_file,
          tumor_ssls_stan_data,
          iter_warmup = iter_warmup,
          iter_sampling = iter_sampling,
          save_warmup = TRUE,
          parallel_chains = chains,
          chains = chains,
          threads_per_chain = tumor_ssls_stan_data$n_shards,
          init = tumor_ssls_initializer,
          adapt_delta = 0.8,
          save_metric = TRUE,
          output_dir = file.path(publication_output_path, "fit", str_c(base_name, "_", disease)),
          timestamp = fit_output_timestamp
        ),
        storage = "main",
        resources = tar_resources(
          crew = tar_resources_crew(controller = "fit")
        )
      ),

      tar_target(
        tumor_ssls_nuts_param,
        bayesplot::nuts_params(tumor_ssls_res)
      ),

      tar_target(
        tumor_ssls_nuts_summary,
        summarize_nuts(tumor_ssls_nuts_param)
      ),

      # Draw extraction ----------------------------------------------------------

      tar_target(
        tumor_ssls_draws_pop,
        select_draws(
          tumor_ssls_res,
          ends_with("_pop"),
          starts_with("pop_"),
          measure_sd_sld,
          matches("_sd_level_"),
          matches("^(time_invariant|time_varying)_coef")
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_draws_patient_params,
        select_draws(
          tumor_ssls_res,
          matches("^(frac|init|tr)_.+_patient"),
          matches("patient_log_(growth|decrease)_rate")
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_draws_sld_recist,
        select_draws(
          tumor_ssls_res,
          rep_patient_log_sld,
          forecast_patient_log_sld,
          rep_recist,
          forecast_obs_recist
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_draws_endpoints,
        select_draws(
          tumor_ssls_res,
          matches("(spop|sample)(_target|_ms)?_(((quant_)?(pfs|os))|km_est|right_censored|(pfs|os)_n)"),
          matches("(spop|sample)_target_(((un)?confirmed_response)|orr)"),
          matches("(spop|sample)_(os|pfs)_(quant|km_est|n)"),
          matches("(spop|sample)_os(_censored)?"),
          matches("(spop|sample)_(os|pfs)_quant_exceeds_max"),
          recist_confusion_matrix
        ),
        resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
      ),

      tar_target(
        tumor_ssls_convergence,
        check_convergence(tumor_ssls_draws_pop, tumor_ssls_draws_patient_params)
      ),

      # KM curves ----------------------------------------------------------------

      tar_target(
        tumor_ssls_km_rvar,
        tumor_ssls_draws_endpoints |>
          recover_types(select(all_analysis_data, trial)) |>
          spread_rvars(
            sample_target_km_est[trial, t],
            spop_target_km_est[trial, t],
            sample_ms_pfs_km_est[trial, t],
            spop_ms_pfs_km_est[trial, t],
            sample_pfs_km_est[trial, t],
            spop_pfs_km_est[trial, t]
          ) |>
          mutate(fit_type = type)
      ),

      tar_target(
        tumor_ssls_km_os_rvar,
        tumor_ssls_draws_endpoints |>
          recover_types(select(all_analysis_data, trial)) |>
          spread_rvars(
            sample_os_km_est[trial, t],
            spop_os_km_est[trial, t]
          ) |>
          mutate(fit_type = type)
      ),

      tar_target(
        tumor_ssls_trial_pfs_quant,
        tumor_ssls_draws_endpoints |>
          recover_types(select(all_analysis_data, trial)) |>
          spread_rvars(
            sample_target_pfs_quant[trial, q],
            spop_target_pfs_quant[trial, q],
            sample_ms_pfs_quant[trial, q],
            spop_ms_pfs_quant[trial, q],
            sample_pfs_quant[trial, q],
            spop_pfs_quant[trial, q]
          ) |>
          left_join(
            enframe(tumor_ssls_stan_data$pfs_quantiles, name = "q", value = "quantile"),
            by = "q"
          ) |>
          mutate(fit_type = type)
      ),

      tar_target(
        tumor_ssls_orr_rvar,
        tumor_ssls_draws_endpoints |>
          recover_types(select(all_analysis_data, trial)) |>
          spread_rvars(
            sample_target_orr[trial],
            spop_target_orr[trial]
          ) |>
          mutate(fit_type = type)
      ),

      tar_target(
        tumor_ssls_forecast_target_pfs_n_rvar,
        tumor_ssls_draws_endpoints |>
          recover_types(select(all_analysis_data, trial)) |>
          spread_rvars(
            sample_target_pfs_n[trial, n],
            spop_target_pfs_n[trial, n],
            sample_ms_pfs_n[trial, n],
            spop_ms_pfs_n[trial, n],
            sample_pfs_n[trial, n],
            spop_pfs_n[trial, n]
          ) |>
          left_join(pfs_timepoints_pub, by = "n") |>
          mutate(fit_type = type)
      ),

      # Population parameters ----------------------------------------------------

      tar_target(
        tumor_ssls_rates_rvar,
        gather_rvars(
          tumor_ssls_draws_pop,
          tr_loc_pop,
          frac_logit_loc_pop,
          tr_sd_level_intercept,
          pop_log_decrease_rate,
          pop_log_growth_rate
        ) |>
          mutate(.value_exp = exp(.value), fit_type = type)
      ),

      tar_target(
        tumor_ssls_noise_sd_rvar,
        gather_rvars(tumor_ssls_draws_pop, measure_sd_sld) |>
          mutate(fit_type = type)
      ),

      tar_map(
        tibble(level = c("patient")),
        names = "level",

        tar_target(
          tumor_ssls_rates_bpi,
          get_tumor_ssls_level_param_binned(
            tumor_ssls_draws_patient_params,
            level,
            param = str_c(
              "{level}_",
              c("log_decrease_rate", "log_growth_rate",
                "log_growth_rate_residual", "log_decrease_rate_residual")
            ),
            type,
            breaks = seq(-3, 3, 0.05),
            inv_link_breaks = seq(0, 25, 0.5)
          )
        ),

        tar_target(
          tumor_ssls_decrease_prop_bpi,
          get_tumor_ssls_level_param_binned(
            tumor_ssls_draws_patient_params,
            level,
            param = str_c("frac_logit_loc_{level}"),
            type,
            breaks = seq(-5, 5, 0.05),
            inv_link = rvar_plogis,
            inv_link_breaks = seq(0, 1, 0.01)
          )
        )
      ),

      tar_map(
        tibble(
          event_type   = c("right_censored", "uncensored"),
          event_cond   = c(expr(right_censored), expr(!right_censored)),
          event_slicer = c(\(d, n) d, \(d, n) slice_sample(d, n = n))
        ),
        names = "event_type",

        tar_target(
          state_patient_subsample,
          get_state_patients(
            all_analysis_data,
            by = trial,
            cond = event_cond,
            slicer = event_slicer,
            sample_size = 40
          )
        ),

        tar_target(
          tumor_ssls_rep_sld_rvar,
          get_sld(tumor_ssls_draws_sld_recist, state_patient_subsample) |>
            mutate(fit_type = type),
          resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
        ),

        tar_target(
          tumor_ssls_forecast_sld_rvar,
          get_forecast_sld(
            tumor_ssls_draws_sld_recist,
            all_analysis_data,
            state_patient_subsample,
            forecast_extent = extend_max_all_t
          ) |>
            mutate(fit_type = type),
          resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
        ),

        tar_target(
          tumor_ssls_staged_sld_rvar,
          bind_rows(
            obs      = tumor_ssls_rep_sld_rvar |> rename(patient_sld = rep_patient_sld),
            forecast = tumor_ssls_forecast_sld_rvar |> rename(patient_sld = forecast_patient_sld),
            .id = "stage"
          )
        ),

        tar_target(
          tumor_ssls_recist_rvar,
          get_recist(tumor_ssls_draws_sld_recist, state_patient_subsample) |>
            mutate(fit_type = type),
          resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
        ),

        tar_target(
          tumor_ssls_forecast_recist_rvar,
          get_forecast_recist(
            tumor_ssls_draws_sld_recist,
            all_analysis_data,
            state_patient_subsample,
            forecast_extent = extend_max_all_t
          ) |>
            mutate(fit_type = type),
          resources = tar_resources(crew = tar_resources_crew(controller = "many samples"))
        ),

        tar_target(
          tumor_ssls_staged_recist_rvar,
          bind_rows(
            obs      = tumor_ssls_recist_rvar |>
              filter(!is.na(response)) |>
              rename(recist = rep_recist),
            forecast = tumor_ssls_forecast_recist_rvar |>
              rename(recist = forecast_obs_recist),
            .id = "stage"
          )
        )
      ),

      tar_target(
        tumor_ssls_coef,
        if (base_tumor_ssls_stan_data$n_covar > 0) {
          gather_rvars(
            tumor_ssls_draws_pop,
            frac_coef_qr_pop[n],
            init_coef_qr_pop[n],
            time_invariant_coef_qr_01[n],
            time_varying_coef_01[n]
          ) |>
            mutate(fit_type = type, .exp_value = exp(.value))
        }
      )
    ),

    # Combined prior + posterior -------------------------------------------------

    tar_target(
      all_tumor_ssls_km_rvar,
      bind_rows(tumor_ssls_km_rvar_prior, tumor_ssls_km_rvar_posterior)
    ),
    tar_target(
      all_tumor_ssls_km_os_rvar,
      bind_rows(tumor_ssls_km_os_rvar_prior, tumor_ssls_km_os_rvar_posterior)
    ),
    tar_target(
      all_tumor_ssls_trial_pfs_quant,
      bind_rows(tumor_ssls_trial_pfs_quant_prior, tumor_ssls_trial_pfs_quant_posterior)
    ),
    tar_target(
      all_tumor_ssls_orr_rvar,
      bind_rows(tumor_ssls_orr_rvar_prior, tumor_ssls_orr_rvar_posterior)
    ),
    tar_target(
      all_tumor_ssls_forecast_target_pfs_n_rvar,
      bind_rows(
        tumor_ssls_forecast_target_pfs_n_rvar_prior,
        tumor_ssls_forecast_target_pfs_n_rvar_posterior
      )
    ),
    tar_target(
      all_tumor_ssls_rates_rvar,
      bind_rows(tumor_ssls_rates_rvar_prior, tumor_ssls_rates_rvar_posterior)
    ),
    tar_target(
      all_tumor_ssls_coef,
      bind_rows(tumor_ssls_coef_prior, tumor_ssls_coef_posterior)
    ),
    tar_target(
      all_tumor_ssls_noise_sd_rvar,
      bind_rows(tumor_ssls_noise_sd_rvar_prior, tumor_ssls_noise_sd_rvar_posterior)
    ),
    tar_target(
      all_tumor_ssls_patient_rates_bpi,
      bind_rows(
        tumor_ssls_rates_bpi_patient_prior,
        tumor_ssls_rates_bpi_patient_posterior
      )
    ),
    tar_target(
      all_tumor_ssls_patient_decrease_prop_bpi,
      bind_rows(
        tumor_ssls_decrease_prop_bpi_patient_prior,
        tumor_ssls_decrease_prop_bpi_patient_posterior
      )
    )
  )
)

publication_targets
