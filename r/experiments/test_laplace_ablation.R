# Phase 1b: ABLATION — which term breaks the tumor Laplace approximation?
# ============================================================================
# Runs the built-in Laplace (mode 1) under a grid of conditions and reports the
# divergence rate + R-hat for each. The full test diverged ~91%; here we peel
# back terms to localize the cause:
#
#   variant          n_latent   censoring (LOD)
#   ----------------------------------------------
#   1lat_noLOD          1        none
#   1lat_LOD            1        yes
#   3lat_noLOD          3        none
#   3lat_LOD            3        yes   (== full test)
#
# SPEED: max_treedepth=7 caps leapfrog steps so a pathological chain diverges
# fast instead of stalling 30 min; 2 chains, 300 warmup / 300 sampling, 8 pts.
# We care about the divergence RATE pattern, not precise posteriors.
# ============================================================================

library(cmdstanr)
library(posterior)
library(dplyr)
library(purrr)

set_cmdstan_path("~/.cmdstan/cmdstan-2.39.0")
set.seed(42)

n_patients <- 8
n_visits_per_patient <- 8

true <- list(tr_loc_pop=-3.0, tr_sd=0.3, frac_logit_pop=0.5, frac_sd=0.3,
             init_logit_pop=0.0, init_sd=0.3, measure_sd=0.15)
# Normalized SLD hovers near 1 (normalized to baseline), so a meaningful LOD is
# near that scale: 0.85 censors ~15% of visits (vs 0.08 which censored none).
lod <- 0.85
log_lod <- log(lod)
visit_weeks <- seq(1, by = 4, length.out = n_visits_per_patient)

z_tr_true   <- rnorm(n_patients)
z_frac_true <- rnorm(n_patients)
z_init_true <- rnorm(n_patients)

# Build one observation vector; `censor=TRUE` encodes below-LOD as 0, else keeps
# the positive simulated value (so the same patients, just no normal_lcdf term).
make_obs <- function(censor) {
  obs <- c(); idx <- c(); pos <- integer(n_patients + 1); pos[1] <- 1L
  for (i in seq_len(n_patients)) {
    tr_loc <- true$tr_loc_pop + true$tr_sd * z_tr_true[i]
    frac_logit <- true$frac_logit_pop + true$frac_sd * z_frac_true[i]
    init_logit <- true$init_logit_pop + true$init_sd * z_init_true[i]
    ldf <- plogis(frac_logit, log.p=TRUE); lgf <- plogis(frac_logit, lower.tail=FALSE, log.p=TRUE)
    dr <- exp(tr_loc+ldf); gr <- exp(tr_loc+lgf)
    ild <- plogis(init_logit, log.p=TRUE); ilg <- plogis(init_logit, lower.tail=FALSE, log.p=TRUE)
    for (v in seq_along(visit_weeks)) {
      dt <- visit_weeks[v]-1
      lp <- matrixStats::logSumExp(c(ild-dr*dt, ilg+gr*dt))
      o <- exp(rnorm(1, lp, true$measure_sd))
      if (censor && o < lod) o <- 0
      obs <- c(obs, o); idx <- c(idx, visit_weeks[v])
    }
    pos[i+1] <- pos[i] + n_visits_per_patient
  }
  list(obs=obs, idx=idx, pos=pos)
}

dat_lod   <- make_obs(TRUE)
dat_nolod <- make_obs(FALSE)
cat(sprintf("LOD variant: %d censored of %d\n",
            sum(dat_lod$obs==0), length(dat_lod$obs)))

model <- cmdstan_model("stan/experiments/laplace_tumor_ablation.stan",
                       include_paths=c("stan","stan/tumor"), quiet=TRUE)

pop <- c("tr_loc_pop","tr_sd","frac_logit_pop","frac_sd","init_logit_pop","init_sd")
init_fn <- function() list(
  tr_loc_pop=rnorm(1,-3,0.1), tr_sd=0.3+abs(rnorm(1,0,0.05)),
  frac_logit_pop=rnorm(1,0.5,0.1), frac_sd=0.3+abs(rnorm(1,0,0.05)),
  init_logit_pop=rnorm(1,0,0.1), init_sd=0.3+abs(rnorm(1,0,0.05)))

run_variant <- function(label, n_latent, dat) {
  d <- list(n_patients=n_patients, n_total_visits=length(dat$obs), n_latent=n_latent,
            normalized_obs=dat$obs, patient_visit_pos=dat$pos, visit_time_idx=dat$idx,
            measure_sd=true$measure_sd, log_lod=log_lod, laplace_mode=1L)
  init_mode <- function() init_fn()
  t0 <- proc.time()
  fit <- tryCatch(model$sample(data=d, init=init_mode, seed=123, chains=2,
                    parallel_chains=2, iter_warmup=300, iter_sampling=300,
                    max_treedepth=7, adapt_delta=0.95, refresh=0,
                    show_messages=FALSE),
                  error=function(e) {cat("  ERROR:", conditionMessage(e), "\n"); NULL})
  elapsed <- (proc.time()-t0)["elapsed"]
  if (is.null(fit)) return(tibble(variant=label, n_latent=n_latent, status="errored",
                                  pct_div=NA, max_rhat=NA, min_ess=NA, secs=round(elapsed,1)))
  diag <- fit$diagnostic_summary(quiet=TRUE)
  s <- fit$summary(variables=pop)
  ndraws <- 300*length(diag$num_divergent)
  tibble(variant=label, n_latent=n_latent, status="ok",
         pct_div=round(100*sum(diag$num_divergent)/ndraws,1),
         max_rhat=round(max(s$rhat,na.rm=TRUE),2),
         min_ess=round(min(s$ess_bulk,na.rm=TRUE)),
         secs=round(elapsed,1))
}

cat("\nRunning ablation grid (each: 2 chains, 300+300, treedepth 7)...\n\n")
grid <- list(
  list("1lat_noLOD", 1, dat_nolod),
  list("1lat_LOD",   1, dat_lod),
  list("3lat_noLOD", 3, dat_nolod),
  list("3lat_LOD",   3, dat_lod)
)
res <- map(grid, \(g) run_variant(g[[1]], g[[2]], g[[3]])) |> list_rbind()

cat("\n========================== ABLATION RESULTS ==========================\n")
print(res)
cat("\nReading: high pct_div / rhat>1.1 = Laplace approx fails for that variant.\n")
cat("If 1lat_noLOD is clean but 3lat_LOD diverges, the extra latents and/or\n")
cat("the LOD censoring are the culprits — compare rows to localize.\n")
