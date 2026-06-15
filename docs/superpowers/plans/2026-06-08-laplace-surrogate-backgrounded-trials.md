# Laplace Surrogate for Backgrounded Trials Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the non-log-concave bi-exponential Laplace path for backgrounded historical trials with a quadratic-in-time surrogate (linear in latent coefficients → Gaussian marginal → Laplace exact), tied to the real population rate parameters through a fixed-anchor mechanistic bridge.

**Architecture:** A backgrounded patient's normalized log-burden is modeled as `β_i0 + β_i1·t + β_i2·t²`. Because the latents `β_i` enter linearly, `laplace_marginal_tol` with solver 1 is exact. `(β_pop, Σ_β)` are NOT free: `β_pop` is the quadratic interpolating the exact bi-exponential log-burden at 3 fixed anchor times evaluated at the population rates; `Σ_β = J·diag(sd²)·Jᵀ` propagates `{tr_sd, frac_sd, init_sd}` through the anchor map's autodiff Jacobian. Work is gated: a standalone Stan+R experiment must reproduce the full-HMC population posterior before any change touches `stan/tumor/sf-ssls-lfo.stan`.

**Tech Stack:** Stan (CmdStan 2.39, `laplace_marginal_tol`), CmdStanR, R (tidyverse, posterior), the existing `stan/modules/` `#include` architecture.

---

## File Structure

**Phase 1 — Standalone gate (cheap, fast iteration, no production code):**
- Create: `stan/experiments/laplace_surrogate_test.stan` — toy model, modes 0 (full-HMC bi-exponential reference) and 1 (surrogate-Laplace). Self-contained `functions{}` block (no module includes yet) so it compiles in seconds and is easy to reason about.
- Create: `r/experiments/test_laplace_surrogate.R` — simulate from the true bi-exponential, fit both modes, apply the population-agreement gate.

**Phase 2 — Production module + integration (ONLY after Phase 1 PASSES):**
- Create: `stan/modules/laplace_surrogate/surrogate.stanfunctions` — the three reusable functions promoted verbatim from the validated experiment: `surrogate_anchor_betas`, `surrogate_bridge_cov`, `surrogate_ll` (functor) + `surrogate_K_fn`.
- Modify: `stan/tumor/sf-ssls-lfo.stan` — `#include` the new module; add the backgrounded-patient `laplace_marginal_tol` term in the `model{}` block; add the anchor-times data input.
- Create: `tests/testthat/test-laplace-surrogate.R` — R-side unit tests for the bridge math (Vandermonde solve, Jacobian propagation) computed independently in R and checked against Stan `expose_functions` output.

**Design boundary:** the bridge math lives entirely in `surrogate.stanfunctions` and is exercised by a toy in Phase 1 and unit tests in Phase 2; `sf-ssls-lfo.stan` only *calls* it for backgrounded patients. One new code path, reachable only by patients we never forecast, rejoining the model at `θ_pop`.

---

## Conventions (read before starting)

- **Native pipe `|>`**, never `%>%`. Tidyverse throughout R.
- **No `install.packages()`**, no hardcoded subject IDs (synthetic only), no `<<-`, no `tar_config_set()`.
- **CmdStan path:** `set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")` — 2.38.0 lacks `laplace_marginal_tol`.
- **Honest gate discipline:** check `pct_divergent`/`rhat`/`ess_bulk` FIRST; never let a low-ESS fit fake a PASS via inflated MCSE.
- **Stan syntax check command:** `~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor <file>`

---

## Task 1: Surrogate Stan model — bridge functions + likelihood functor

**Files:**
- Create: `stan/experiments/laplace_surrogate_test.stan`

The whole bridge + surrogate lives in the `functions{}` block of this experiment first, so we validate the exact code we will later promote to the module. Three pieces: (a) `surrogate_anchor_betas` maps population rates → quadratic coefficients via the constant `V⁻¹`; (b) `surrogate_bridge_cov` propagates the SDs through that map's Jacobian; (c) `surrogate_ll` is the `laplace_marginal_tol` functor.

- [ ] **Step 1: Write the model file**

```stan
// Phase 1: log-concave QUADRATIC surrogate for backgrounded-trial patients,
// validated against the full bi-exponential via the built-in laplace_marginal_tol.
// ============================================================================
// Backgrounded patients' normalized log-burden is modeled as a quadratic in t:
//     log_burden_i(t) = beta_i0 + beta_i1*t + beta_i2*t^2 + N(0, measure_sd^2)
// Latents beta_i enter LINEARLY => Gaussian marginal => Laplace is EXACT
// => solver 1 (PD-Hessian Cholesky) is valid (unlike the bi-exponential, which
//    needed solver 3 just to survive).
//
// (beta_pop, Sigma_beta) are NOT free: they are the mechanistic-bridge image of
// the population rate params (tr_loc_pop, frac_logit_pop, init_logit_pop) and
// the per-patient SDs (tr_sd, frac_sd, init_sd):
//   beta_pop = Vinv * [g(t0), g(t1), g(t2)]   (g = exact bi-exponential log-burden)
//   Sigma_beta = J * diag(sd^2) * J'          (J = d beta_pop / d theta_pop, autodiff)
//
// Modes (laplace_mode): 0 = full-HMC bi-exponential reference, 1 = surrogate-Laplace.
// ============================================================================

functions {
  // Exact bi-exponential predicted normalized log-burden at time offset t (weeks).
  // Same rate construction as the full model (sf.stanfunctions / tumor module).
  real bi_exp_log_burden(real t,
                         real tr_loc, real frac_logit, real init_logit) {
    real log_dec_frac = log_inv_logit(frac_logit);
    real log_gro_frac = log1m_inv_logit(frac_logit);
    real dec_rate = exp(tr_loc + log_dec_frac);
    real gro_rate = exp(tr_loc + log_gro_frac);
    real init_log_dec = log_inv_logit(init_logit);
    real init_log_gro = log1m_inv_logit(init_logit);
    real state_dec = init_log_dec - dec_rate * t;
    real state_gro = fmin(init_log_gro + gro_rate * t, 500.0);
    return log_sum_exp(state_dec, state_gro);
  }

  // Quadratic coefficients (beta0, beta1, beta2) that interpolate the exact
  // bi-exponential log-burden at the 3 fixed anchor times. Vinv is the constant
  // inverse Vandermonde passed in as data.
  vector surrogate_anchor_betas(real tr_loc_pop, real frac_logit_pop,
                                real init_logit_pop,
                                data matrix Vinv, data vector anchor_times) {
    vector[3] g;
    for (k in 1:3)
      g[k] = bi_exp_log_burden(anchor_times[k], tr_loc_pop,
                               frac_logit_pop, init_logit_pop);
    return Vinv * g;
  }

  // Analytic Jacobian of the exact bi-exponential log-burden g(t) w.r.t. the
  // population location params (tr_loc, frac_logit, init_logit), at a single t.
  // g = log_sum_exp(a, b), a = init_log_dec - dec_rate*t, b = init_log_gro + gro_rate*t.
  // dg = w_dec*da + w_gro*db with w_dec = softmax weight on the decay branch.
  // (Stan has NO callable autodiff Jacobian of a user function — the `jacobian`
  //  block/`jacobian +=` are for custom transforms, not a returnable matrix — so
  //  we differentiate the closed form directly. It is elementary here.)
  row_vector bi_exp_log_burden_grad(real t, real tr_loc, real frac_logit,
                                    real init_logit) {
    real p = inv_logit(frac_logit);          // dec fraction
    real q = inv_logit(init_logit);           // initial dec share
    real log_dec_frac = log_inv_logit(frac_logit);
    real log_gro_frac = log1m_inv_logit(frac_logit);
    real dec_rate = exp(tr_loc + log_dec_frac);
    real gro_rate = exp(tr_loc + log_gro_frac);
    real a = log_inv_logit(init_logit) - dec_rate * t;
    real b = fmin(log1m_inv_logit(init_logit) + gro_rate * t, 500.0);
    real m = fmax(a, b);
    real w_dec = exp(a - m) / (exp(a - m) + exp(b - m));
    real w_gro = 1 - w_dec;
    // partials of a, b w.r.t. (tr_loc, frac_logit, init_logit)
    real da_dtr  = -t * dec_rate;             real db_dtr  =  t * gro_rate;
    real da_dfr  = -t * dec_rate * (1 - p);   real db_dfr  = -t * gro_rate * p;
    real da_din  = 1 - q;                     real db_din  = -q;
    row_vector[3] g_grad;
    g_grad[1] = w_dec * da_dtr + w_gro * db_dtr;
    g_grad[2] = w_dec * da_dfr + w_gro * db_dfr;
    g_grad[3] = w_dec * da_din + w_gro * db_din;
    return g_grad;
  }

  // Bridge covariance Sigma_beta = J diag(sd^2) J', J = d beta_pop / d theta_pop.
  // beta_pop = Vinv * g(anchors), so J = Vinv * [grad g(t_k)]_k (a 3x3 stack of
  // the per-anchor gradients). A diagonal jitter keeps K PD if an sd -> 0.
  matrix surrogate_bridge_cov(real tr_loc_pop, real frac_logit_pop,
                              real init_logit_pop,
                              real tr_sd, real frac_sd, real init_sd,
                              data matrix Vinv, data vector anchor_times,
                              data real jitter) {
    matrix[3, 3] g_jac;   // row k = d g(t_k) / d theta
    for (k in 1:3)
      g_jac[k] = bi_exp_log_burden_grad(anchor_times[k], tr_loc_pop,
                                        frac_logit_pop, init_logit_pop);
    matrix[3, 3] J = Vinv * g_jac;
    matrix[3, 3] D = diag_matrix(square([tr_sd, frac_sd, init_sd]'));
    return J * D * J' + diag_matrix(rep_vector(jitter, 3));
  }

  // laplace_marginal_tol functor. theta = stacked per-patient [b0,b1,b2] latents
  // in the WHITENED coordinate (K = identity, mean 0); we map to the bridged
  // distribution inside via the Cholesky of Sigma_beta passed through phi.
  // Here we use the simpler route: K_fn returns Sigma_beta directly, so theta is
  // in the natural beta coordinate centered at beta_pop. So the functor receives
  // beta_pop and adds (theta_i - 0) usage: latents are deviations from beta_pop.
  real surrogate_ll(vector theta,
                    vector beta_pop,
                    real measure_sd, real log_lod,
                    data int n_patients,
                    data vector normalized_obs,
                    data array[] int patient_visit_pos,
                    data array[] int visit_time,
                    data array[] int patient_of_visit) {
    int d = 3;
    real lp = 0;
    for (i in 1:n_patients) {
      real b0 = beta_pop[1] + theta[(i - 1) * d + 1];
      real b1 = beta_pop[2] + theta[(i - 1) * d + 2];
      real b2 = beta_pop[3] + theta[(i - 1) * d + 3];
      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real t = visit_time[v] - 1.0;
        real mu = b0 + b1 * t + b2 * t * t;
        if (normalized_obs[v] > 0)
          lp += normal_lpdf(log(normalized_obs[v]) | mu, measure_sd);
        else
          lp += normal_lcdf(log_lod | mu, measure_sd);
      }
    }
    return lp;
  }

  // Prior covariance functor for laplace_marginal_tol: the bridged Sigma_beta,
  // block-replicated per patient is handled by hessian_block_size=3 + this K
  // returning the per-block covariance tiled. We return a full block-diagonal.
  matrix surrogate_K_fn(matrix Sigma_beta, int n_patients) {
    int d = 3;
    matrix[n_patients * d, n_patients * d] K =
      rep_matrix(0, n_patients * d, n_patients * d);
    for (i in 1:n_patients) {
      int s = (i - 1) * d + 1;
      K[s:(s + d - 1), s:(s + d - 1)] = Sigma_beta;
    }
    return K;
  }
}

data {
  int<lower=1> n_patients;
  int<lower=1> n_total_visits;
  vector[n_total_visits] normalized_obs;       // <= 0 marks below-LOD
  array[n_patients + 1] int patient_visit_pos;
  array[n_total_visits] int visit_time;        // 1-indexed week (t = idx - 1)
  array[n_total_visits] int patient_of_visit;  // patient id per visit row
  real<lower=0> measure_sd;
  real log_lod;
  vector[3] anchor_times;                      // fixed calendar anchors (weeks)
  int<lower=0, upper=1> laplace_mode;
}

transformed data {
  int latent_dim = n_patients * 3;
  vector[latent_dim] theta_0 = rep_vector(0.0, latent_dim);
  real tolerance = 1e-8;
  int max_num_steps = 100;        // log-concave => Newton converges fast
  int hessian_block_size = 3;
  int solver = 1;                 // PD-Hessian Cholesky: valid for log-concave
  int max_steps_line_search = 0;  // not needed when well-conditioned
  int allow_fallback = 1;
  real jitter = 1e-8;

  // Constant inverse Vandermonde for the 3 fixed anchors.
  matrix[3, 3] V;
  for (k in 1:3) {
    V[k, 1] = 1.0;
    V[k, 2] = anchor_times[k];
    V[k, 3] = anchor_times[k] * anchor_times[k];
  }
  matrix[3, 3] Vinv = inverse(V);
}

parameters {
  real tr_loc_pop;
  real<lower=0> tr_sd;
  real frac_logit_pop;
  real<lower=0> frac_sd;
  real init_logit_pop;
  real<lower=0> init_sd;

  vector[laplace_mode == 0 ? n_patients : 0] z_tr;
  vector[laplace_mode == 0 ? n_patients : 0] z_frac;
  vector[laplace_mode == 0 ? n_patients : 0] z_init;
}

model {
  tr_loc_pop     ~ normal(-3, 1);
  tr_sd          ~ normal(0, 1);
  frac_logit_pop ~ normal(0.5, 1);
  frac_sd        ~ normal(0, 1);
  init_logit_pop ~ normal(0, 1);
  init_sd        ~ normal(0, 1);

  if (laplace_mode == 0) {
    // MODE 0: full-HMC bi-exponential reference (the true generative model)
    z_tr ~ std_normal();
    z_frac ~ std_normal();
    z_init ~ std_normal();
    for (i in 1:n_patients) {
      real tr_loc     = tr_loc_pop     + tr_sd   * z_tr[i];
      real frac_logit = frac_logit_pop + frac_sd * z_frac[i];
      real init_logit = init_logit_pop + init_sd * z_init[i];
      int v_start = patient_visit_pos[i];
      int v_end   = patient_visit_pos[i + 1] - 1;
      for (v in v_start:v_end) {
        real t = visit_time[v] - 1.0;
        real mu = bi_exp_log_burden(t, tr_loc, frac_logit, init_logit);
        if (normalized_obs[v] > 0)
          target += normal_lpdf(log(normalized_obs[v]) | mu, measure_sd);
        else
          target += normal_lcdf(log_lod | mu, measure_sd);
      }
    }
  } else {
    // MODE 1: surrogate-Laplace (the thing under test)
    vector[3] beta_pop = surrogate_anchor_betas(
      tr_loc_pop, frac_logit_pop, init_logit_pop, Vinv, anchor_times);
    matrix[3, 3] Sigma_beta = surrogate_bridge_cov(
      tr_loc_pop, frac_logit_pop, init_logit_pop,
      tr_sd, frac_sd, init_sd, Vinv, anchor_times, jitter);

    target += laplace_marginal_tol(
      surrogate_ll,
      (beta_pop, measure_sd, log_lod, n_patients,
       normalized_obs, patient_visit_pos, visit_time, patient_of_visit),
      hessian_block_size,
      surrogate_K_fn,
      (Sigma_beta, n_patients),
      (theta_0, tolerance, max_num_steps, solver, max_steps_line_search, allow_fallback)
    );
  }
}
```

- [ ] **Step 2: Syntax-check the model**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/experiments/laplace_surrogate_test.stan
```
Expected: no output (clean parse).

**Note on the Jacobian:** the spec said "built-in autodiff Jacobian," but Stan has **no callable autodiff Jacobian of a user function** — the `jacobian` block and `jacobian +=` statement are for declaring custom-transform log-det adjustments, not for returning `∂f/∂x` as a matrix inside a function. So `bi_exp_log_burden_grad` differentiates the closed form analytically (it's elementary: `g = log_sum_exp(a,b)` ⇒ `∂g = w_dec·∂a + w_gro·∂b` with softmax weights). The result is identical to what autodiff would produce, just hand-derived. **Validation safety net:** Task 4's unit test finite-differences `surrogate_anchor_betas` in R and checks the analytic `J` against it, so a derivative typo cannot pass silently.

- [ ] **Step 3: Commit**

```bash
git add stan/experiments/laplace_surrogate_test.stan
git commit -m "feat(laplace): quadratic surrogate Stan model for backgrounded trials (gate)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: R validation harness — simulate, fit both modes, apply gate

**Files:**
- Create: `r/experiments/test_laplace_surrogate.R`

Simulate backgrounded patients from the TRUE bi-exponential (with LOD censoring), fit mode 0 (reference) and mode 1 (surrogate), and gate on population-posterior agreement.

- [ ] **Step 1: Write the harness**

```r
# Phase 1 gate: does the quadratic surrogate (Laplace, solver 1) reproduce the
# population posterior of the full bi-exponential (full HMC) for backgrounded
# patients? Data is simulated from the TRUE bi-exponential, so the surrogate is
# tested against the real generative process, not its own assumptions.
# ============================================================================
library(cmdstanr)
library(posterior)
library(dplyr)
library(tibble)
library(purrr)

set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(42)

n_patients <- 12
n_visits_per_patient <- 8

true <- list(tr_loc_pop = -3.0, tr_sd = 0.3, frac_logit_pop = 0.5,
             frac_sd = 0.3, init_logit_pop = 0.0, init_sd = 0.3,
             measure_sd = 0.15)
lod <- 0.85               # ~15% censored at the normalized-to-baseline scale
log_lod <- log(lod)
visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)

# Fixed calendar anchors (weeks): baseline, mid, late — span the visit window.
anchor_times <- c(0, 12, 28)

z_tr   <- rnorm(n_patients)
z_frac <- rnorm(n_patients)
z_init <- rnorm(n_patients)

obs <- c(); idx <- c(); pid <- c()
pos <- integer(n_patients + 1); pos[1] <- 1L
for (i in seq_len(n_patients)) {
  tr_loc     <- true$tr_loc_pop     + true$tr_sd   * z_tr[i]
  frac_logit <- true$frac_logit_pop + true$frac_sd * z_frac[i]
  init_logit <- true$init_logit_pop + true$init_sd * z_init[i]
  ldf <- plogis(frac_logit, log.p = TRUE)
  lgf <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
  dr <- exp(tr_loc + ldf); gr <- exp(tr_loc + lgf)
  ild <- plogis(init_logit, log.p = TRUE)
  ilg <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)
  for (v in seq_along(visit_weeks)) {
    t <- visit_weeks[v] - 1
    lp <- matrixStats::logSumExp(c(ild - dr * t, ilg + gr * t))
    o <- exp(rnorm(1, lp, true$measure_sd))
    if (o < lod) o <- 0
    obs <- c(obs, o); idx <- c(idx, visit_weeks[v]); pid <- c(pid, i)
  }
  pos[i + 1] <- pos[i] + n_visits_per_patient
}

cat(sprintf("Data: %d patients, %d visits, %d censored (%.0f%%)\n",
            n_patients, length(obs), sum(obs == 0),
            100 * sum(obs == 0) / length(obs)))

stan_data <- list(
  n_patients = n_patients, n_total_visits = length(obs),
  normalized_obs = obs, patient_visit_pos = pos, visit_time = idx,
  patient_of_visit = pid, measure_sd = true$measure_sd, log_lod = log_lod,
  anchor_times = anchor_times)

model <- cmdstan_model("stan/experiments/laplace_surrogate_test.stan",
                       include_paths = c("stan", "stan/tumor"), quiet = FALSE)

pop <- c("tr_loc_pop", "tr_sd", "frac_logit_pop", "frac_sd",
         "init_logit_pop", "init_sd")
init_fn <- function() list(
  tr_loc_pop = rnorm(1, -3, 0.1), tr_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  frac_logit_pop = rnorm(1, 0.5, 0.1), frac_sd = 0.3 + abs(rnorm(1, 0, 0.05)),
  init_logit_pop = rnorm(1, 0, 0.1), init_sd = 0.3 + abs(rnorm(1, 0, 0.05)))

fits <- list(); timings <- list()
for (mode in 0:1) {
  name <- c("full_hmc", "surrogate")[mode + 1]
  cat(sprintf("\n=== mode %d (%s) ===\n", mode, name))
  init_mode <- if (mode == 0) {
    function() c(init_fn(), list(z_tr = rnorm(n_patients, 0, 0.5),
                                 z_frac = rnorm(n_patients, 0, 0.5),
                                 z_init = rnorm(n_patients, 0, 0.5)))
  } else init_fn
  t0 <- proc.time()
  fits[[name]] <- model$sample(
    data = c(stan_data, list(laplace_mode = mode)), init = init_mode,
    seed = 123, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 2000, refresh = 500,
    adapt_delta = if (mode == 1) 0.95 else 0.9)
  timings[[name]] <- (proc.time() - t0)["elapsed"]
}

# --- Population comparison -------------------------------------------------
hmc <- fits$full_hmc$summary(variables = pop)
sur <- fits$surrogate$summary(variables = pop)
comparison <- tibble(
  parameter = hmc$variable, hmc_mean = hmc$mean, surrogate_mean = sur$mean,
  diff = abs(hmc$mean - sur$mean), hmc_sd = hmc$sd, surrogate_sd = sur$sd,
  mcse = hmc$sd / sqrt(pmin(hmc$ess_bulk, sur$ess_bulk)),
  diff_over_mcse = abs(hmc$mean - sur$mean) /
    (hmc$sd / sqrt(pmin(hmc$ess_bulk, sur$ess_bulk))),
  sd_ratio = sur$sd / hmc$sd)
cat("\n=== population comparison ===\n"); print(comparison)
cat(sprintf("\nTimings: HMC=%.1fs  surrogate=%.1fs\n",
            timings$full_hmc, timings$surrogate))

# --- Joint posterior: correlation-structure check --------------------------
# Subsumes target-forecast invariance (forecast is a nonlinear fn of the joint).
cor_hmc <- fits$full_hmc$draws(variables = pop, format = "draws_matrix") |>
  cor()
cor_sur <- fits$surrogate$draws(variables = pop, format = "draws_matrix") |>
  cor()
max_cor_diff <- max(abs(cor_hmc - cor_sur))
cat(sprintf("\nMax |corr difference| (joint structure): %.3f\n", max_cor_diff))

# --- Convergence gate (FIRST) ----------------------------------------------
diag_one <- function(fit, name) {
  d <- fit$diagnostic_summary(quiet = TRUE)
  s <- fit$summary(variables = pop)
  n_draws <- fit$metadata()$iter_sampling * length(d$num_divergent)
  tibble(mode = name, pct_divergent = 100 * sum(d$num_divergent) / n_draws,
         max_rhat = max(s$rhat, na.rm = TRUE),
         min_ess_bulk = min(s$ess_bulk, na.rm = TRUE))
}
diag_tbl <- bind_rows(diag_one(fits$full_hmc, "full_hmc"),
                      diag_one(fits$surrogate, "surrogate"))
cat("\n=== convergence diagnostics ===\n"); print(diag_tbl)

sur_d <- diag_tbl |> filter(mode == "surrogate")
converged <- sur_d$pct_divergent < 1 && sur_d$max_rhat < 1.01 &&
  sur_d$min_ess_bulk > 400
max_ratio <- max(comparison$diff_over_mcse)
min_sd_ratio <- min(comparison$sd_ratio)

cat(sprintf("\nMax diff/MCSE: %.2f | min SD ratio (sur/hmc): %.2f | max corr diff: %.3f\n",
            max_ratio, min_sd_ratio, max_cor_diff))
if (!converged) {
  cat("NO-GO: surrogate fit did not converge (solver 1 should be clean if",
      "the surrogate is log-concave as designed).\n")
} else if (max_ratio < 5 && min_sd_ratio > 0.8 && max_cor_diff < 0.1) {
  cat("PASS: surrogate converged AND population posterior (means, SDs, joint",
      "correlations) agrees with full HMC.\n")
} else {
  cat("FAIL: surrogate converged but population posterior differs",
      "(approximation bias) — see which of diff/MCSE, SD ratio, corr diff broke.\n")
}
cat("\nPhase 1 gate complete.\n")
```

- [ ] **Step 2: Run the gate**

Run:
```bash
Rscript r/experiments/test_laplace_surrogate.R
```
Expected: prints data summary, both fits complete, a `PASS` / `FAIL` / `NO-GO` verdict. The decisive lines are the convergence table (surrogate `pct_divergent < 1`, `max_rhat < 1.01`, `min_ess_bulk > 400`) and then the agreement metrics.

- [ ] **Step 3: Commit the harness and record the verdict**

```bash
git add r/experiments/test_laplace_surrogate.R
git commit -m "test(laplace): population-agreement gate for quadratic surrogate

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## GATE — STOP HERE if Phase 1 did not PASS

**Do not start Task 3 until Task 2 prints `PASS`.**

- **PASS** → proceed to Task 3 (promote to module + integrate).
- **FAIL (converged but biased)** → the quadratic is too stiff for this data's curvature. Apply spec §6 fallback: extend the surrogate to a cubic / 4-knot natural spline (still linear-in-coefficients → Laplace stays exact). Concretely: widen `anchor_times` to 4 knots, change `[3]`/`3` dims to `[4]`/`4` throughout Task 1's file, rebuild the `V`/`Vinv` as 4×4, extend `hessian_block_size` and the latent stride to 4. Re-run Task 2. Update the spec's risk section with the outcome.
- **NO-GO (surrogate did not converge under solver 1)** → contradicts the core design claim (linear-Gaussian ⇒ log-concave). Stop and re-derive; do not integrate. Record findings in `laplace-builtin-2-39.md`.

Record the verdict and key numbers (divergences, diff/MCSE, SD ratio, corr diff, timings) in the memory file before continuing.

---

## Task 3: Promote validated functions to a reusable module

**Files:**
- Create: `stan/modules/laplace_surrogate/surrogate.stanfunctions`

Move the four validated functions (`bi_exp_log_burden`, `surrogate_anchor_betas`, `surrogate_bridge_cov`, `surrogate_ll`, `surrogate_K_fn`) out of the experiment file into a module include, unchanged except for doc comments. This is the DRY step: one definition, used by both the experiment (optional) and the production model.

- [ ] **Step 1: Create the module file** by copying the exact, validated bodies of the five functions from `stan/experiments/laplace_surrogate_test.stan`'s `functions{}` block into `stan/modules/laplace_surrogate/surrogate.stanfunctions`. Prepend a module header docstring describing the bridge (mirror spec §2–§3). Do NOT alter the function bodies — they are the validated artifact.

- [ ] **Step 2: Syntax-check the module via a throwaway includer**

Run:
```bash
printf 'functions {\n#include "modules/laplace_surrogate/surrogate.stanfunctions"\n}\n' > /tmp/_inc_check.stan
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor /tmp/_inc_check.stan
```
Expected: clean parse (no output).

- [ ] **Step 3: Commit**

```bash
git add stan/modules/laplace_surrogate/surrogate.stanfunctions
git commit -m "feat(laplace): extract validated surrogate bridge into reusable module

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: R unit tests for the bridge math

**Files:**
- Create: `tests/testthat/test-laplace-surrogate.R`

Independently compute `β_pop` and `Σ_β` in R and check Stan's exposed functions match. This guards the bridge math against silent regressions when the module is later touched.

- [ ] **Step 1: Write the test**

```r
# Unit tests for the surrogate bridge math (Vandermonde solve + Jacobian cov).
# We expose the Stan functions and compare against an independent R computation.
library(testthat)
library(cmdstanr)
set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")

test_that("surrogate_anchor_betas reproduces an R Vandermonde interpolation", {
  # Tiny includer model that exposes just the bridge functions.
  stan_file <- write_stan_file(paste0(
    "functions {\n",
    "#include \"modules/laplace_surrogate/surrogate.stanfunctions\"\n",
    "}\nmodel {}\n"))
  m <- cmdstan_model(stan_file, include_paths = c("stan", "stan/tumor"),
                     compile_standalone = TRUE)
  m$expose_functions(global = TRUE)

  anchors <- c(0, 12, 28)
  tr_loc <- -3.0; frac_logit <- 0.5; init_logit <- 0.0
  V <- cbind(1, anchors, anchors^2)
  Vinv <- solve(V)

  bi_exp_r <- function(t) {
    ldf <- plogis(frac_logit, log.p = TRUE)
    lgf <- plogis(frac_logit, lower.tail = FALSE, log.p = TRUE)
    dr <- exp(tr_loc + ldf); gr <- exp(tr_loc + lgf)
    ild <- plogis(init_logit, log.p = TRUE)
    ilg <- plogis(init_logit, lower.tail = FALSE, log.p = TRUE)
    matrixStats::logSumExp(c(ild - dr * t, ilg + gr * t))
  }
  g <- vapply(anchors, bi_exp_r, numeric(1))
  beta_expected <- as.numeric(Vinv %*% g)

  beta_stan <- surrogate_anchor_betas(tr_loc, frac_logit, init_logit,
                                      Vinv, anchors)
  expect_equal(beta_stan, beta_expected, tolerance = 1e-8)
})

test_that("analytic bridge Jacobian matches a finite-difference of beta_pop", {
  # Guards bi_exp_log_burden_grad against a hand-derivation typo.
  anchors <- c(0, 12, 28); V <- cbind(1, anchors, anchors^2); Vinv <- solve(V)
  theta0 <- c(-3.0, 0.5, 0.0)
  beta_at <- function(th) surrogate_anchor_betas(th[1], th[2], th[3], Vinv, anchors)
  eps <- 1e-5
  J_fd <- vapply(seq_along(theta0), function(j) {
    tp <- theta0; tp[j] <- tp[j] + eps
    tm <- theta0; tm[j] <- tm[j] - eps
    (beta_at(tp) - beta_at(tm)) / (2 * eps)
  }, numeric(3))                                   # 3x3: rows=beta, cols=theta
  # Reconstruct analytic J = Vinv %*% g_jac via the cov function's internals:
  # Sigma = J diag(sd^2) J'; with sd=1 and no jitter, Sigma = J J'. Instead we
  # expose J directly by setting sd=(1,0,0) etc. to read columns — simpler: trust
  # surrogate_bridge_cov with sd=1,jitter=0 gives J J', and compare to J_fd J_fd'.
  Sigma_analytic <- surrogate_bridge_cov(theta0[1], theta0[2], theta0[3],
                                         1, 1, 1, Vinv, anchors, 0)
  expect_equal(Sigma_analytic, J_fd %*% diag(3) %*% t(J_fd), tolerance = 1e-4)
})

test_that("surrogate_bridge_cov is symmetric positive definite", {
  anchors <- c(0, 12, 28); V <- cbind(1, anchors, anchors^2); Vinv <- solve(V)
  Sigma <- surrogate_bridge_cov(-3.0, 0.5, 0.0, 0.3, 0.3, 0.3,
                                Vinv, anchors, 1e-8)
  expect_equal(Sigma, t(Sigma), tolerance = 1e-10)   # symmetric
  expect_true(all(eigen(Sigma, only.values = TRUE)$values > 0))  # PD
})
```

- [ ] **Step 2: Run the test**

Run:
```bash
Rscript -e 'testthat::test_file("tests/testthat/test-laplace-surrogate.R")'
```
Expected: all three tests PASS. The bridge functions are plain Stan (no `jacobian`/no higher-order autodiff), so `expose_functions` compiles them cleanly. `surrogate_ll`/`surrogate_K_fn` need not be exposed — they're exercised end-to-end by Task 2's gate.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-laplace-surrogate.R
git commit -m "test(laplace): unit tests for surrogate bridge (Vandermonde + cov PD)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Integrate the surrogate into `sf-ssls-lfo.stan`

**Files:**
- Modify: `stan/tumor/sf-ssls-lfo.stan` — functions include (after line 12), a data input for anchor times, and the backgrounded-patient term in `model{}` (after line 100).

Backgrounded patients currently contribute nothing to the LFO likelihood (their NCP params don't exist — params are sized by `n_forecast_patients`; see `stan/_full_model_transformed_data.stan:100-114`). We add their surrogate-marginalized contribution. Forecast patients (the `normalized_sld ~ sf_log_space_obs` term at line 98) are untouched.

- [ ] **Step 1: Add the module include**

In `stan/tumor/sf-ssls-lfo.stan`, after line 12 (`#include "modules/tumor/tumor.stanfunctions"`), add:

```stan
  #include "modules/laplace_surrogate/surrogate.stanfunctions"
```

- [ ] **Step 2: Add the anchor-times data input**

In the `data {}` block (after line 35, `#include "modules/state_space/lfo_data.stan"`), add:

```stan
  // Fixed calendar anchors (weeks) for the backgrounded-trial Laplace surrogate.
  // Set in R to span the backgrounded trials' observation window (spec §6).
  vector[3] surrogate_anchor_times;
  int<lower=0, upper=1> enable_background_surrogate;
```

- [ ] **Step 3: Build the constant inverse Vandermonde in transformed data**

In the `transformed data {}` block (after line 57, `#include "_lfo_transformed_data.stan"`), add:

```stan
  // Constant inverse Vandermonde + Laplace control constants for the surrogate.
  matrix[3, 3] surrogate_V;
  for (k in 1:3) {
    surrogate_V[k, 1] = 1.0;
    surrogate_V[k, 2] = surrogate_anchor_times[k];
    surrogate_V[k, 3] = surrogate_anchor_times[k] * surrogate_anchor_times[k];
  }
  matrix[3, 3] surrogate_Vinv = inverse(surrogate_V);
  int surrogate_latent_dim = n_background_patients * 3;
  vector[surrogate_latent_dim] surrogate_theta_0 =
    rep_vector(0.0, surrogate_latent_dim);
  // Visit rows and positions for backgrounded patients only, in a compact array.
  // (Reuse patient_visit_pos with background_patient_idx; build a local compact
  //  index so the functor sees a contiguous 1..n_background_patients layout.)
```

- [ ] **Step 4: Add the backgrounded-patient likelihood term**

In the `model {}` block, inside `if (fit_tumor_data) { ... }`, after the multistate block (after line 109's closing `}` for `enable_ms_01`, before line 110's closing `}`), add:

```stan
    // --- Backgrounded historical trials: surrogate-marginalized SLD ---
    // These patients inform the population rate params only; their per-patient
    // latents are integrated out via the log-concave quadratic surrogate.
    if (enable_background_surrogate && n_background_patients > 0) {
      // Build compact (1..n_background_patients) views of the obs/positions.
      array[n_background_patients + 1] int bg_pos;
      bg_pos[1] = 1;
      for (j in 1:n_background_patients) {
        int p = background_patient_idx[j];
        int vs, ve;
        (vs, ve) = get_pos(patient_visit_pos, p);
        bg_pos[j + 1] = bg_pos[j] + (ve - vs + 1);
      }
      int n_bg_visits = bg_pos[n_background_patients + 1] - 1;
      vector[n_bg_visits] bg_obs;
      array[n_bg_visits] int bg_time;
      array[n_bg_visits] int bg_patient;
      {
        int w = 1;
        for (j in 1:n_background_patients) {
          int p = background_patient_idx[j];
          int vs, ve;
          (vs, ve) = get_pos(patient_visit_pos, p);
          for (v in vs:ve) {
            // normalized to baseline; below-LOD encoded as <= 0 as in functor
            bg_obs[w] = normalized_sld[v];
            bg_time[w] = t_patient_visit_idx[v];
            bg_patient[w] = j;
            w += 1;
          }
        }
      }

      vector[3] bg_beta_pop = surrogate_anchor_betas(
        tr_loc_pop, frac_logit_pop, init_logit_pop,
        surrogate_Vinv, surrogate_anchor_times);
      matrix[3, 3] bg_Sigma = surrogate_bridge_cov(
        tr_loc_pop, frac_logit_pop, init_logit_pop,
        tr_sd, frac_sd, init_sd, surrogate_Vinv, surrogate_anchor_times, 1e-8);

      target += laplace_marginal_tol(
        surrogate_ll,
        (bg_beta_pop, measure_sd_sld, log_lod, n_background_patients,
         bg_obs, bg_pos, bg_time, bg_patient),
        3,                            // hessian_block_size
        surrogate_K_fn,
        (bg_Sigma, n_background_patients),
        (surrogate_theta_0, 1e-8, 100, 1, 0, 1)  // tol, max_steps, solver=1, ls, fallback
      );
    }
```

**Note on parameter names:** `tr_loc_pop`, `frac_logit_pop`, `init_logit_pop`, `tr_sd`, `frac_sd`, `init_sd`, `measure_sd_sld`, `log_lod` must match the names in scope inside this model's `model{}` block (they come from the `tr`/`frac`/`init`/`state_space` modules). Before writing, grep to confirm the exact in-scope symbol names and fix any mismatch:
```bash
grep -rn "tr_loc_pop\|frac_logit_pop\|init_logit_pop\|measure_sd_sld\|real log_lod" stan/modules/tr stan/modules/frac stan/modules/init stan/modules/state_space
```
If a population symbol is actually an array/hierarchy element (e.g. `tr_loc_pop` is a level-0 intercept), use the population-level accessor that the existing forecast path uses for the same quantity — match it exactly.

- [ ] **Step 5: Syntax-check the full model**

Run:
```bash
~/.cmdstan/cmdstan-2.39.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssls-lfo.stan
```
Expected: clean parse. Fix any name/scope errors surfaced (this is where Step 4's name-check pays off).

- [ ] **Step 6: Commit**

```bash
git add stan/tumor/sf-ssls-lfo.stan
git commit -m "feat(laplace): marginalize backgrounded-trial latents via quadratic surrogate

Forecast patients unchanged; backgrounded patients contribute a
solver-1 Laplace-marginalized SLD term tied to population rates through
the fixed-anchor bridge. Gated behind enable_background_surrogate.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Wire the new data inputs into the R pipeline

**Files:**
- Modify: `r/priors.R` — add `surrogate_anchor_times` default (and `enable_background_surrogate` default off).
- Modify: the tumor stan-data assembly (find with grep below) — pass the two new fields.

- [ ] **Step 1: Locate the stan-data assembly and prior defaults**

Run:
```bash
grep -rn "log_lod\|fit_tumor_data\|n_forecast_patients" r/ --include=*.R | grep -iv test | head
grep -n "measure_sd_sld_alpha\|tr_loc_pop_mean" r/priors.R
```
This identifies the function that builds the tumor Stan data list and the prior-list constructor.

- [ ] **Step 2: Add defaults in `r/priors.R`**

In the returned prior/data list, add (placed near the other tumor SLD entries):

```r
    # Backgrounded-trial Laplace surrogate (off by default; anchors span the
    # historical trials' observation window in weeks — tune per run, spec §6)
    enable_background_surrogate = 0L,
    surrogate_anchor_times = c(0, 12, 28),
```

- [ ] **Step 3: Ensure the fields reach the Stan data list**

If the stan-data assembly spreads the prior list (common pattern), Step 2 is sufficient. Otherwise add the two fields explicitly in the data-list constructor identified in Step 1. Confirm with:
```bash
grep -rn "surrogate_anchor_times\|enable_background_surrogate" r/
```
Expected: appears in both the prior defaults and (transitively or explicitly) the stan-data list.

- [ ] **Step 4: Syntax-check the R**

Run:
```bash
Rscript -e 'invisible(parse("r/priors.R"))'
```
Expected: no error.

- [ ] **Step 5: Commit**

```bash
git add r/priors.R
git commit -m "feat(laplace): default surrogate anchors + off-by-default flag in priors

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: End-to-end smoke + memory update

**Files:**
- Modify: `/home/ubuntu/.claude/projects/-mnt-code/memory/laplace-builtin-2-39.md`

- [ ] **Step 1: Compile the production model through CmdStanR** (catches include/data issues the bare `stanc` check misses)

Run:
```bash
Rscript -e 'library(cmdstanr); set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0"); cmdstan_model("stan/tumor/sf-ssls-lfo.stan", include_paths=c("stan","stan/tumor"), compile=TRUE)'
```
Expected: compiles without error.

- [ ] **Step 2: Update the memory file** — flip the NO-GO verdict to the surrogate resolution: record that the bi-exponential Laplace is unusable but the quadratic surrogate (solver 1) reproduces the population posterior (cite the Task 2 gate numbers), and that it's wired behind `enable_background_surrogate` in `sf-ssls-lfo.stan`. Keep the historical NO-GO context.

- [ ] **Step 3: Commit** (memory dir is outside the repo — commit only repo files; the memory file is saved in place by the Write tool).

```bash
git add docs/superpowers/plans/2026-06-08-laplace-surrogate-backgrounded-trials.md
git commit -m "docs(laplace): implementation plan for backgrounded-trial surrogate

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Done criteria

- Phase 1 gate prints `PASS` (surrogate clean under solver 1; population means within 5 MCSE, SD ratio > 0.8, max corr diff < 0.1).
- Module functions extracted and unit-tested.
- `sf-ssls-lfo.stan` compiles with the surrogate term, gated behind `enable_background_surrogate` (default off — no behavior change until a run opts in and sets anchors).
- Memory file reflects the surrogate resolution.
