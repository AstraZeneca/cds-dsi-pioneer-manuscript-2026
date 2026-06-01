export const meta = {
  name: 'review-spop-pfs-spec',
  description: 'Focused independent review of the spop dropout PFS-from-OS design spec across correctness / CIF bit-identity / statistical-convention lenses, synthesized into a go/no-go verdict',
  phases: [
    { title: 'Critique', detail: '3 independent critics, one per risk lens' },
    { title: 'Synthesize', detail: 'merge findings into go/no-go + concrete revisions' },
  ],
}

const WT = (args && args.worktree) || '/mnt/code/.worktrees/karim/msm-frailty-impl'
const SPEC = `${WT}/docs/superpowers/specs/2026-05-30-spop-dropout-pfs-convention-design.md`

const base = `You are reviewing a design spec at ${SPEC}. READ IT FIRST (it is short). ` +
  `You may also read the referenced code in ${WT}/stan/pfs.stanfunctions, ` +
  `${WT}/stan/multistate.stanfunctions, ${WT}/stan/tumor/_tumor_endpoints_generated_quantities.stan, ` +
  `and ${WT}/r/util.R to ground your critique. Read-only — do NOT edit anything. ` +
  `The spec changes how spop PFS is derived for 0->3 dropout patients: PFS becomes the OS ` +
  `outcome (event at 3->2 death, else censored at horizon), via reordering the OS draw before ` +
  `the PFS derivation and grafting (approach A), with the CIF cause-attribution re-keyed off the ` +
  `explicit spop_is_dropout flag.`

const FINDING = {
  type: 'object',
  properties: {
    lens: { type: 'string' },
    verdict: { type: 'string', enum: ['approve', 'approve-with-changes', 'block'] },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          severity: { type: 'string', enum: ['blocker', 'major', 'minor', 'nit'] },
          issue: { type: 'string' },
          where: { type: 'string', description: 'spec section or file:line' },
          suggested_fix: { type: 'string' },
        },
        required: ['severity', 'issue', 'suggested_fix'],
      },
    },
  },
  required: ['lens', 'verdict', 'findings'],
}

// ── Phase 1: three independent critics, distinct lenses (barrier: synthesis needs all) ──
phase('Critique')
const critics = await parallel([
  () => agent(`${base}\n\nLENS: STAN GQ CORRECTNESS. Verify the approach-A reorder is sound: ` +
    `(a) does derive_spop_os_rng truly have no dependency on PFS outputs, so moving it earlier is safe? ` +
    `(b) does grafting spop_os into spop_pfs for cause 3 preserve PFS<=OS by construction? ` +
    `(c) are there cause-3 sub-cases (e.g. SLD target PD before dropout, the old case at pfs.stanfunctions:1395) ` +
    `that the blanket graft would wrongly override? (d) sample-path untouched — confirmed correct? ` +
    `Cite file:line. Default to flagging anything uncertain.`,
    { label: 'critic:stan-correctness', phase: 'Critique', schema: FINDING }),
  () => agent(`${base}\n\nLENS: CIF BIT-IDENTITY. The spec re-keys the CIF drop_/dd/prog predicates off ` +
    `spop_is_dropout in BOTH the Stan GQ and the r/util.R mirror. Verify: (a) is spop_is_dropout actually ` +
    `available at both sites with the same semantics? (b) will the proposed predicate rewrite produce ` +
    `identical Stan and R results across all cause combinations (progression, on-trial death, dropout-death, ` +
    `dropout-alive, fully censored)? (c) any draw where the old and new attribution differ for a NON-dropout ` +
    `patient (which would be an unintended change)? Cite file:line.`,
    { label: 'critic:cif-identity', phase: 'Critique', schema: FINDING }),
  () => agent(`${base}\n\nLENS: STATISTICAL / CLINICAL CONVENTION. Does PFS-from-OS for dropouts actually ` +
    `match the observed-data convention (pfs = death_week-1 for death-without-progression; censor at last ` +
    `assessment otherwise)? Specifically: (a) the observed cohort censors a non-dying dropout at last visit ` +
    `(~wk21) but the spec censors at max_all_t — is comparing these on one KM legitimate, or does it bias the ` +
    `spop optimistic? (b) is treating an off-trial death as a PFS *event* defensible under standard PFS ` +
    `definitions, or should it be censored? (c) does the change risk double-counting deaths between the PFS ` +
    `and OS curves? Reason as a trial statistician.`,
    { label: 'critic:stat-convention', phase: 'Critique', schema: FINDING }),
])
const found = critics.filter(Boolean)

// ── Phase 2: synthesize into a single go/no-go ──
phase('Synthesize')
const synthesis = await agent(
  `${base}\n\nYou are the lead reviewer. Three independent critics returned these structured findings:\n` +
  JSON.stringify(found, null, 2) +
  `\n\nSynthesize into a single verdict. Deduplicate overlapping findings. Rank by severity. ` +
  `Give an overall go/no-go: APPROVE (spec is sound as-is), REVISE (list the specific spec edits needed ` +
  `before planning), or BLOCK (a fundamental flaw — explain). Be concrete: each required revision should ` +
  `name the spec section and the exact change. Do NOT edit the spec — just report.`,
  { label: 'synthesis', phase: 'Synthesize', schema: {
      type: 'object',
      properties: {
        overall: { type: 'string', enum: ['approve', 'revise', 'block'] },
        summary: { type: 'string' },
        required_revisions: {
          type: 'array',
          items: {
            type: 'object',
            properties: {
              section: { type: 'string' },
              change: { type: 'string' },
              rationale: { type: 'string' },
            },
            required: ['section', 'change'],
          },
        },
        nonblocking_notes: { type: 'array', items: { type: 'string' } },
      },
      required: ['overall', 'summary', 'required_revisions'],
    } })

return { critics: found, synthesis }
