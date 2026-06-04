export const meta = {
  name: 'implement-correlated-msm-frailty',
  description: 'Implement decomposed MSM level effects + correlated frailty intercepts per the 2026-05-29 design spec, gated through degeneracy + prior-predictive checks',
  whenToUse: 'After the correlated-MSM-frailty design spec is approved; stops before the multi-hour Domino fit so it can be launched interactively.',
  phases: [
    { title: 'Discovery', detail: 'parallel read-only mapping of every refactor touch-point' },
    { title: 'Refactor', detail: 'serial: decompose shared flag into 3 per-transition arrays (no correlation yet)' },
    { title: 'Gate1', detail: 'adversarial degeneracy + syntax verification of the refactor' },
    { title: 'Frailty', detail: 'serial: add LKJ-Cholesky correlated intercept block' },
    { title: 'Gate2', detail: 'syntax + dimensioning + validation + NCP-assembly review' },
    { title: 'PriorPredictive', detail: 'configure {01,03} patient correlation, MCMC-free spop CIF plausibility check' },
  ],
}

// ---------------------------------------------------------------------------
// Shared context. WORKTREE is created by the operator BEFORE launch (see the
// plan message) and passed via args.worktree. Every editing/gate agent works
// in that ONE shared tree so the dependent serial chain stays coherent — we do
// NOT use per-agent isolation:'worktree' (that would give each agent its own
// disconnected copy and break Phase 2 building on Phase 1).
// ---------------------------------------------------------------------------
const WT = (args && args.worktree) || '/mnt/code/.worktrees/karim/pioneer-pub'
const SPEC = `${WT}/docs/superpowers/specs/2026-05-29-correlated-msm-frailty-design.md`
const MOD = `${WT}/stan/modules/multistate`

const ctx = `You are implementing the design at ${SPEC}. READ IT FIRST.
Work ONLY in the worktree ${WT}. The multistate module is at ${MOD}.
Key facts from the design:
- Replace the single shared flag enable_ms_level_baseline_hazard[lv] (0-4) with three
  per-transition x per-level arrays: enable_ms_level_gp, ms_level_intercept_mode,
  ms_level_intercept_corr_group.
- Intercept slots: {01,02,03,12_s,12_t,32}. Only {01,03} are correlated for this analysis.
- Correlated block: cholesky_factor_corr L + NCP std-normal z; u = diag_pre_multiply(sigma,L)*z.
  sigma reuses existing per-transition level intercept SDs. LKJ(2). Gaussian marginals only.
- Singletons / corr_group=0 MUST keep the exact current scalar sigma*raw path => bit-identical.
- NO change to the GQ simulator (calculate_all_patients_endpoints_rng) or r/util.R.
- Follow stan-guidelines.md (NCP, zeros_vector, fatal_error for validation) and
  the 8-step "Adding Module Parameters" checklist.`

const FINDINGS = {
  type: 'object',
  properties: {
    sites: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          lines: { type: 'string' },
          what: { type: 'string', description: 'what must change here' },
        },
        required: ['file', 'lines', 'what'],
      },
    },
    notes: { type: 'string' },
  },
  required: ['sites'],
}

const VERDICT = {
  type: 'object',
  properties: {
    pass: { type: 'boolean' },
    evidence: { type: 'string', description: 'command output or diff excerpt proving the verdict' },
    issues: { type: 'array', items: { type: 'string' } },
  },
  required: ['pass', 'evidence'],
}

// ===========================================================================
// PHASE 0 — Discovery (parallel fan-out, read-only). Genuine barrier: the
// serial edit phase needs the COMPLETE map before it starts.
// ===========================================================================
phase('Discovery')
const discovery = await parallel([
  () => agent(`${ctx}\n\nMAP TASK: trace enable_ms_level_baseline_hazard through ${MOD}/transformed_data.stan and the group-count helpers it calls (compute_ms_transition_group_counts, split_cp_ncp_pos, compute_ms_level_baseline_flags). List every sizing expression that must generalize from one shared mode vector to per-transition mode vectors. Read-only.`,
    { label: 'map:sizing', phase: 'Discovery', schema: FINDINGS }),
  () => agent(`${ctx}\n\nMAP TASK: inventory all raw/cp_log_lambda_gp_0k_level_intercept declarations and GP eta blocks in ${MOD}/parameters.stan and ${MOD}/hyperparams.stan, for all 6 intercept slots {01,02,03,12_s,12_t,32}. Read-only.`,
    { label: 'map:params', phase: 'Discovery', schema: FINDINGS }),
  () => agent(`${ctx}\n\nMAP TASK: locate every "log_cond_surv_0k += rep_matrix(...)" intercept-addition site in ${MOD}/transformed_parameters.stan (expect one per slot). These are where the frailty rides in. Read-only.`,
    { label: 'map:tp', phase: 'Discovery', schema: FINDINGS }),
  () => agent(`${ctx}\n\nMAP TASK: find the flag and its consumers in ${WT}/r/priors.R, ${WT}/r/initializers.R, ${WT}/targets/publication_targets.R, and ${WT}/tests/testthat/. Identify where the config-translation helper (legacy flag -> 3 arrays) should live. Read-only.`,
    { label: 'map:r', phase: 'Discovery', schema: FINDINGS }),
])
const workList = discovery.filter(Boolean).flatMap(d => d.sites)
log(`Discovery mapped ${workList.length} touch-points across Stan + R.`)

// ===========================================================================
// PHASE 1 — Decomposition refactor (SERIAL, single agent). One coherent diff;
// shared group-count contracts make this unsplittable. NO correlation yet.
// ===========================================================================
phase('Refactor')
const refactor = await agent(
  `${ctx}\n\nIMPLEMENT PHASE 1 ONLY (decomposition, NO correlated block yet).\n` +
  `Complete work-list from discovery:\n${JSON.stringify(workList, null, 2)}\n\n` +
  `Steps: (1) add the 3 new arrays to flags.stan; (2) generalize transformed_data group-sizing ` +
  `per-transition; (3) split GP-residual sizing from intercept sizing; (4) update parameters/` +
  `transformed_parameters/priors/hyperparams so each transition reads its own mode; (5) add the ` +
  `R config-translation helper so the legacy single-flag value reproduces bit-identically; ` +
  `(6) wire publication_targets.R to pass the 3 arrays (still all-independent, no corr_group>0). ` +
  `Run stanc syntax check on stan/tumor/sf-ssm-log-space.stan before returning. ` +
  `Report exactly which files changed and the stanc result.`,
  { label: 'refactor:decompose', phase: 'Refactor' })
log('Refactor phase complete. Verifying degeneracy...')

// ----- Gate 1: adversarial degeneracy verification (parallel skeptics) -----
phase('Gate1')
const gate1 = await parallel([
  () => agent(`${ctx}\n\nVERIFY: run stanc (with --include-paths=stan --include-paths=stan/tumor) on stan/tumor/sf-ssm-log-space.stan AND stanc (--include-paths=stan --include-paths=stan/psa) on stan/psa/pioneer.stan in ${WT}. pass=true ONLY if BOTH parse clean. Paste stanc output as evidence.`,
    { label: 'gate1:syntax', phase: 'Gate1', schema: VERDICT }),
  () => agent(`${ctx}\n\nVERIFY: run the testthat suite focused on multistate dimensioning/degeneracy in ${WT} (Rscript -e 'testthat::test_dir("tests/testthat")' or the targeted dimensioning test file). pass=true ONLY if the no-correlation config produces output identical to the legacy single-flag model. Paste test summary as evidence.`,
    { label: 'gate1:tests', phase: 'Gate1', schema: VERDICT }),
  () => agent(`${ctx}\n\nADVERSARIAL REVIEW: your explicit job is to REFUTE the claim "the decomposition is bit-identical when no correlation is configured." Read the refactor diff (git diff in ${WT}). Hunt for any place where per-transition sizing diverges from the legacy shared-flag sizing, any group-count off-by-one, any intercept dropped or double-counted. Default to pass=false if you find ANY plausible behavior change. Cite file:line as evidence.`,
    { label: 'gate1:skeptic', phase: 'Gate1', schema: VERDICT }),
])
const g1 = gate1.filter(Boolean)
const g1pass = g1.length === 3 && g1.every(v => v.pass)
if (!g1pass) {
  log('GATE 1 FAILED — refactor changed behavior. Issues: ' +
    JSON.stringify(g1.flatMap(v => v.issues || [])))
  return { stoppedAt: 'Gate1', refactor, gate1: g1 }
}
log('GATE 1 PASSED — refactor is degeneracy-safe.')

// ===========================================================================
// PHASE 2 — Correlated block (SERIAL, single agent). Small, additive.
// ===========================================================================
phase('Frailty')
const frailty = await agent(
  `${ctx}\n\nIMPLEMENT PHASE 2 (correlated intercept block) on top of the verified refactor.\n` +
  `Add per-(level,group) cholesky_factor_corr L_ms_intercept_corr + NCP matrix z_ms_intercept to ` +
  `parameters.stan; MVN assembly u=diag_pre_multiply(sigma_members,L)*z in transformed_parameters.stan ` +
  `scattered to each member transition's level intercept via patient_ms_baseline_flat_idx; ` +
  `lkj_corr_cholesky(2) + std_normal() priors in priors.stan; fatal_error validation in ` +
  `transformed_data.stan (corr_group member must be RE; members share group partition; singleton collapses). ` +
  `Add LKJ eta to r/priors.R get_multistate_priors() and identity-Cholesky + rnorm(sd=0.3) inits to ` +
  `r/initializers.R. Keep singletons on the exact scalar path. stanc-check before returning.`,
  { label: 'frailty:add-block', phase: 'Frailty' })

// ----- Gate 2: correctness (parallel) -----
phase('Gate2')
const gate2 = await parallel([
  () => agent(`${ctx}\n\nVERIFY: stanc on sf-ssm-log-space.stan and pioneer.stan in ${WT}. pass only if both clean. Paste output.`,
    { label: 'gate2:syntax', phase: 'Gate2', schema: VERDICT }),
  () => agent(`${ctx}\n\nVERIFY: run testthat in ${WT}; the new dimensioning + validation tests for the correlated block must pass and the all-zero-corr degeneracy test must STILL pass. Paste summary.`,
    { label: 'gate2:tests', phase: 'Gate2', schema: VERDICT }),
  () => agent(`${ctx}\n\nADVERSARIAL REVIEW: scrutinize the NCP assembly u=diag_pre_multiply(sigma,L)*z — correct dimension order, sigma_members really the per-transition level SDs, z is std_normal, scatter maps row m to the right transition's intercept vector. Confirm singletons (corr_group=0) bypass the Cholesky entirely. Default pass=false on any doubt. Cite file:line.`,
    { label: 'gate2:ncp-review', phase: 'Gate2', schema: VERDICT }),
])
const g2 = gate2.filter(Boolean)
const g2pass = g2.length === 3 && g2.every(v => v.pass)
if (!g2pass) {
  log('GATE 2 FAILED. Issues: ' + JSON.stringify(g2.flatMap(v => v.issues || [])))
  return { stoppedAt: 'Gate2', frailty, gate2: g2 }
}
log('GATE 2 PASSED — correlated block correct.')

// ===========================================================================
// PHASE 3 — Prior-predictive gate (single agent, MCMC-free). Terminal phase.
// ===========================================================================
phase('PriorPredictive')
const priorPred = await agent(
  `${ctx}\n\nFINAL PHASE: configure {01,03} patient-level correlation in ` +
  `targets/publication_targets.R (intercept_mode RE on 01 and 03 at the patient level, ` +
  `corr_group code shared, 02 independent). Then run an MCMC-FREE prior-predictive check: ` +
  `sample L (LKJ2), z (std_normal), sigma from their priors and push through to the spop CIF / ` +
  `dropout-vs-progression routing for the 497 patients. Confirm the frailty produces a PLAUSIBLE ` +
  `spread of 0->3 vs 0->1 routing (not everyone dropping out, not no one). Report the routing summary ` +
  `and whether it passes the plausibility bar. Do NOT launch any Domino job.`,
  { label: 'prior-predictive', phase: 'PriorPredictive', schema: {
      type: 'object',
      properties: {
        plausible: { type: 'boolean' },
        routing_summary: { type: 'string' },
        next_step: { type: 'string', description: 'what the operator should do to launch the fit' },
      },
      required: ['plausible', 'routing_summary'],
    } })

log(priorPred && priorPred.plausible
  ? 'PRIOR-PREDICTIVE PASSED — ready for interactive fit launch.'
  : 'PRIOR-PREDICTIVE flagged implausible routing — review before fitting.')

return {
  worktree: WT,
  refactor, frailty, priorPred,
  stoppedAt: 'PriorPredictive (by design)',
  handoff: 'Launch the Domino posterior fit + downstream (-D) interactively via pioneer-toolkit:start-job, then run $diagnostic_summary() and the died_off_trial spop-vs-observed KM comparison.',
}
