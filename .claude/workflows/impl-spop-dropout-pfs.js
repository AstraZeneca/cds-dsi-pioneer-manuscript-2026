export const meta = {
  name: 'impl-spop-dropout-pfs',
  description: 'Implement the spop dropout PFS-from-OS convention per the 2026-05-30 plan: serial Stan GQ edits (tasks 2-6) then an adversarial verification gate. Stops before the Domino downstream regeneration.',
  phases: [
    { title: 'Implement', detail: 'serial Stan edits across pfs.stanfunctions + endpoints GQ (plan tasks 2-6, 9)' },
    { title: 'Gate', detail: 'adversarial verification: stanc + degeneracy + plan-fidelity review' },
  ],
}

const WT = (args && args.worktree) || '/mnt/code/.worktrees/karim/msm-frailty-impl'
const PLAN = `${WT}/docs/superpowers/plans/2026-05-30-spop-dropout-pfs-convention.md`
const SPEC = `${WT}/docs/superpowers/specs/2026-05-30-spop-dropout-pfs-convention-design.md`

const ctx = `You are implementing the plan at ${PLAN} (spec: ${SPEC}). READ THE PLAN FIRST. ` +
  `Work ONLY in the worktree ${WT}. This is a generated-quantities + R/docs change — NO MCMC re-fit. ` +
  `Key facts:\n` +
  `- The change makes a spop 0->3 dropout's PFS equal its OS outcome (event at 3->2 death, else ` +
  `censored at horizon), by reordering the OS draw before the PFS derivation in the GQ loop and ` +
  `grafting spop_os into spop_pfs/spop_ms_pfs for spop_cause==3.\n` +
  `- spop_dropout_week is threaded through calculate_all_patients_endpoints_rng's return tuple so the ` +
  `3->2 sojourn KM can be recovered after the graft makes spop_pfs==spop_os.\n` +
  `- compute_trial_cif must test is_dropout FIRST so grafted dropout-deaths stay in CIF_03.\n` +
  `- BOTH sojourn blocks need fixing: spop_km_32 (use spop_dropout_week) and spop_km_12 (exclude dropouts).\n` +
  `- The R mirror compute_cif_from_draws was already deleted; do NOT look for it.\n` +
  `- Sample/conditional path (derive_sample_*) is UNTOUCHED.\n` +
  `Stan syntax check command: ~/.cmdstan/cmdstan-2.38.0/bin/stanc --include-paths=stan --include-paths=stan/tumor stan/tumor/sf-ssm-log-space.stan`

// ===========================================================================
// PHASE 1 — Serial implementation (tasks 2-6 share pfs.stanfunctions; task 9 docs).
// ONE agent does the whole edit sequence; the tuple-threading + graft + CIF
// reorder are interdependent and must be internally consistent.
// ===========================================================================
phase('Implement')
const impl = await agent(
  `${ctx}\n\nIMPLEMENT plan Tasks 2 through 6 (the Stan changes) and Task 9 (the docs note), ` +
  `in order, committing after each task exactly as the plan specifies. Do NOT do Task 1 (baseline ` +
  `read), Task 7 (test file — a separate agent owns verification), or Task 8 (Domino regeneration — ` +
  `out of scope for this workflow). \n\n` +
  `CRITICAL correctness points to get exactly right:\n` +
  `1. Task 2: the new spop_dropout_week entry must be the LAST element in BOTH the return-TYPE tuple ` +
  `and the return-VALUE tuple of calculate_all_patients_endpoints_rng — positions must match or Stan ` +
  `mis-binds silently. Verify by counting elements in both.\n` +
  `2. Task 3: the OS-draw block must move ABOVE the derive_spop_pfs call; the graft block goes AFTER ` +
  `derive_spop_pfs; spop_dropout_week[j]=spop_time_03 is captured for cause 3 before the graft.\n` +
  `3. Task 6: the call-site LHS destructuring tuple must gain spop_dropout_week as its final element, ` +
  `matching Task 2's return position.\n` +
  `Run the stanc syntax check after EACH Stan task and confirm clean parse before committing. ` +
  `Report every file changed, each commit SHA, and the final stanc result.`,
  { label: 'implement:stan+docs', phase: 'Implement' })

log('Implementation phase complete. Verifying...')

// ===========================================================================
// PHASE 2 — Adversarial verification gate (parallel, independent).
// ===========================================================================
phase('Gate')
const VERDICT = {
  type: 'object',
  properties: {
    pass: { type: 'boolean' },
    evidence: { type: 'string', description: 'command output or diff excerpt proving the verdict' },
    issues: { type: 'array', items: { type: 'string' } },
  },
  required: ['pass', 'evidence'],
}

const gate = await parallel([
  () => agent(`${ctx}\n\nVERIFY (syntax): run stanc on stan/tumor/sf-ssm-log-space.stan AND ` +
    `stan/psa/pioneer.stan (--include-paths=stan --include-paths=stan/psa for the latter). ` +
    `pass=true ONLY if BOTH parse clean. Also write the Task 7 test file ` +
    `tests/testthat/test-spop-dropout-pfs.R EXACTLY as specified in the plan and commit it (the impl ` +
    `agent skipped it). Paste stanc output as evidence.`,
    { label: 'gate:syntax+test', phase: 'Gate', schema: VERDICT }),
  () => agent(`${ctx}\n\nADVERSARIAL REVIEW — your job is to REFUTE that the implementation matches ` +
    `the plan. Read the git diff in ${WT} for the Stan files. Check specifically: (a) is ` +
    `spop_dropout_week the final element in BOTH the return-type and return-value tuples of ` +
    `calculate_all_patients_endpoints_rng, AND the final element of the call-site LHS destructuring in ` +
    `_tumor_endpoints_generated_quantities.stan? A position mismatch is a silent SEVERE bug. ` +
    `(b) does the OS draw now precede derive_spop_pfs, with the graft AFTER it? (c) does compute_trial_cif ` +
    `test is_dropout FIRST? (d) is spop_km_32 using spop_dropout_week (not spop_pfs), and does spop_km_12 ` +
    `exclude dropouts in BOTH its filter and its n_prog count? (e) is the sample/conditional path truly ` +
    `untouched? Default to pass=false on ANY doubt. Cite file:line.`,
    { label: 'gate:plan-fidelity', phase: 'Gate', schema: VERDICT }),
])
const g = gate.filter(Boolean)
const passed = g.length === 2 && g.every(v => v.pass)
log(passed ? 'GATE PASSED — ready for downstream regeneration.'
           : 'GATE FAILED: ' + JSON.stringify(g.flatMap(v => v.issues || [])))

return {
  implementation: impl,
  gate: g,
  passed,
  handoff: passed
    ? 'Stan + docs + test complete and verified. NEXT (interactive): rebuild downstream targets ' +
      'from the existing fit (generate_quantities / -D -m tumor_ssls_draws_endpoints_posterior, no re-fit), ' +
      'then run plan Task 8 validation: died_off_trial spop median PFS toward 42, CIF_03 count == n_dropout, ' +
      'PFS<=OS, spop_km_32 not collapsed, and the test passes against regenerated draws.'
    : 'Gate failed — see issues; fix and re-run before downstream regeneration.',
}
