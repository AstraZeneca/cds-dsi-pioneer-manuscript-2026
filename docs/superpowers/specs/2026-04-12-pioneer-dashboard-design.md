# Pioneer Dashboard — Design Specification

## Overview

A Shiny-based web dashboard deployed as a Domino App for monitoring pioneer pipeline runs. Combines automated Stan MCMC diagnostics with AI-powered interpretation (Haiku), providing a persistent history of fit quality across commits and model variants.

**Scope (MVP):** Diagnostics + live job status. No job launching from the UI.

**Audience:** Single user (Karim).

**Project:** Pioneer only (extensible to sclc later).

## Architecture

Two-process deployment behind a single Domino App entry point.

```
run.sh (Domino App entry point)
├── Rscript worker.R &       ← background, loops every 15 min
└── Rscript -e 'shiny::runApp(port=8888)'  ← foreground, always on
```

### Worker (worker.R)

Background process with two polling loops:

**Diagnostic loop (every 15 min):**
1. Discover fit directories under `/mnt/data/.../pioneer/*/fit/`
2. Call `stan-diagnose-helper.sh` per target (cache-aware, skips unchanged CSVs via `.cache` mtime comparison)
3. Write numeric diagnostics to SQLite
4. Call Haiku API with numbers + recent history to generate plain-English narrative
5. Store narrative in SQLite

**Job poller (every 2 min):**
1. Read `job-log.jsonl` for new entries
2. Call Domino jobs API for live status of running jobs
3. Write to SQLite `job_events` / update Job records

### Shiny App (app.R)

Reads SQLite only — never touches Stan CSV files. Auto-refreshes via `reactiveTimer` every 60 seconds.

**Panels:**
- Jobs panel: recent launches, live status badges, log tail on click
- Runs panel: stores x model variants, build time, last updated
- Diagnostics panel: Rhat/ESS heatmap per target, divergences, E-BFMI
- History panel: trend charts over time, Haiku narratives timeline

### Shared State

SQLite database at `/mnt/data/analysis-results/karim_naguib/pioneer-dashboard/db.sqlite` — persists across Domino App restarts.

## Domain Model

The domain model combines three design patterns:

1. **Composite** — tree structure mirroring nested `tar_map()` calls
2. **Type Object / Meta-Object Protocol (MOP)** — Variant instances act as types for their children, recursively
3. **Stratified DAG** — typed build dependencies with pipeline stage constraints

### Entity Hierarchy

```
Target (abstract, «component»)
│   id · parent_id · kind · commit_id
│
├── Variant (concrete, «composite · type-object»)
│     store_id · dimension · value · config
│
└── LeafTarget (abstract)
      job_id · target_name · built_at
      │
      ├── UpstreamTarget (concrete)
      │     (no extra fields)
      │     Examples: stan_data, initializer, exe_hash
      │
      ├── FitTarget (concrete)
      │     fit_dir_path · csv_prefix · csv_run_id · model · type
      │     → has many DiagnosticSnapshots
      │
      └── DownstreamTarget (concrete)
            (no extra fields)
            Examples: km_rvar, draws_pop, cond_pfs
```

### Supporting Entities

```
Store
  id · name · path · created_at

Commit
  id · sha · branch · message · authored_at · author

Job
  id · store_id (FK→Store) · commit_id (FK→Commit)
  domino_job_id · targets_requested · launched_at · completed_at · status · hw_tier

DiagnosticSnapshot
  id · fit_target_id (FK→FitTarget)
  captured_at · csv_mtime · rhat_max · ess_bulk_min · n_divergent_total · severity

ChainDiagnostic
  id · snapshot_id (FK→DiagnosticSnapshot)
  chain_num · n_divergent · stepsize · e_bfmi · mean_accept · n_max_treedepth

Narrative
  id · snapshot_id (FK→DiagnosticSnapshot)
  text · created_at · prompt_tokens · severity
```

### Relationships

| From | To | Cardinality | Notes |
|------|-----|-------------|-------|
| Commit | Job | 1:∞ | A commit can be used by multiple jobs |
| Commit | Target | 1:∞ | On base class. Variant = config version; Leaf = build version |
| Store | Job | 1:∞ | Jobs run sequentially against a store |
| Store | Variant | 1:∞ | Only Variants need direct store link |
| Job | LeafTarget | 1:∞ | Each leaf was built by exactly one job |
| Variant | Target (children) | 1:0..* | Self-referential via parent_id. Composite containment + type-instance |
| FitTarget | DiagnosticSnapshot | 1:∞ | Diagnostic history over time |
| DiagnosticSnapshot | ChainDiagnostic | 1:4 | One per MCMC chain |
| DiagnosticSnapshot | Narrative | 1:0..1 | Haiku may not have run yet |

**Constraint:** `leaf.commit_id == leaf.job.commit_id` (denormalized for query speed)

**Leaf targets derive Store through:** LeafTarget → Job → Store (no direct FK needed)

### Composite + Type Object (MOP)

The Variant class serves dual roles simultaneously:

1. **Composite node** — structural container in the tree (holds children via `parent_id`)
2. **Type Object** — carries configuration (`config` field: JSON of tribble columns) that characterizes all its children

When a Variant holds another Variant as a child, the relationship is:
- **Structural containment** (Composite pattern)
- **Type-instance** (child is an instance of the parent's type)
- **Metatype** (the child is itself a type for its own children)

**Meta-Object Levels:**

| Level | Role | Example |
|-------|------|---------|
| M2 (metatype) | Outer tar_map dimension | `Variant(dimension=model, value=combined, config={hist:T, propensity:T, ...})` |
| M1 (type) | Inner tar_map dimension | `Variant(dimension=type, value=posterior)` |
| M0 (instance) | Built artifact | `FitTarget(pioneer_fit_posterior_combined)` |

**Full type of a leaf target** = accumulated config from the entire ancestor Variant chain (M2.config ∪ M1.config ∪ ...).

**Config versioning:** `Variant.commit_id` tracks when the variant's config last changed. When comparing diagnostics across commits, a config change means the targets are potentially "a different kind of thing" — not just rebuilt, but differently configured.

### Stratified Dependency DAG

Build dependencies form a directed acyclic graph on LeafTarget, constrained by pipeline stage:

| Source stage | Valid dependency targets | Constraint |
|-------------|------------------------|------------|
| Upstream | Upstream only | Cannot depend on Fit or Downstream |
| FitTarget | Upstream only | Cannot depend on Fit or Downstream |
| Downstream | Fit + Upstream + Downstream | **Must have ≥1 FitTarget dependency** |

Edges flow strictly forward: Upstream → Fit → Downstream, never backward.

**Staleness detection:** Walk upstream from any target — if any dependency has a newer `commit_id`, this target is stale. Propagates in one direction only.

**Cross-variant dependencies:** A target in one Variant subtree can depend on a target in a different subtree (e.g., shared data preparation targets). The DAG crosses the Composite tree.

### Naming Convention (tar_map)

Target names encode dimensions innermost-first:

```
{base_name}_{type}_{model}
```

Example: `pioneer_fit_posterior_combined_markov_fe_psa_re_ms`

**Outer tar_map (model):** ~20 variants including `combined`, `combined_ungated`, `no_hist`, `arm_only`, `flatiron_only`, `combined_markov_fe_psa_re_ms`, etc.

**Inner tar_map (type):** `posterior`, `prior`

## Worker Details

### stan-diagnose-helper.sh Integration

The worker reuses the existing diagnostic helper script from the bayesian-toolkit plugin. For each fit target:

```bash
bash stan-diagnose-helper.sh <output_dir> <prefix> <run_id> summary
```

The helper:
- Dynamically counts header lines (handles different Stan versions)
- Tail-reads CSVs (avoids loading 10-15GB files into memory)
- Caches results keyed on file modification time
- Returns structured per-chain NUTS statistics

### Haiku AI Layer

After numerical diagnostics are computed, the worker calls Claude Haiku with:

**Input:**
- Current snapshot numbers (Rhat, ESS, divergences, E-BFMI per chain)
- Last N snapshots for this fit target (trend)
- Variant config (what model settings are active)
- Whether config changed since last snapshot

**Output:**
- Plain-English narrative (2-3 sentences)
- Severity classification: `ok` / `warn` / `bad`

**Prompt structure:** Ask Haiku to interpret the diagnostics in context: "Given this is a hierarchical model with [config details], and divergences increased from X to Y since last run at commit Z, what does this suggest?"

### Dependency Graph Population

The worker populates the dependency DAG from `targets::tar_network()` or `tar_manifest()` on each diagnostic loop, updating edges when the pipeline code changes across commits.

## Deployment

### Domino App Configuration

- **Entry point:** `run.sh`
- **Port:** 8888 (Shiny default for Domino Apps)
- **Hardware tier:** Small CPU (diagnostics are read-only, no Stan sampling)
- **Environment:** Must include R, Shiny, RSQLite, httr2 (for Claude API), cmdstanr (for stan-diagnose-helper)

### Data Persistence

- SQLite at `/mnt/data/analysis-results/karim_naguib/pioneer-dashboard/db.sqlite`
- Survives app restarts (on persistent Domino dataset mount)
- Worker creates schema on first run if db doesn't exist

### API Keys

- `ANTHROPIC_API_KEY` or `CLAUDE_API_KEY` — for Haiku calls from the worker
- `DOMINO_USER_API_KEY` — already available in Domino env, for jobs API polling

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Web framework | R Shiny |
| Background worker | R (standalone Rscript process) |
| Database | SQLite via RSQLite / DBI |
| Diagnostics engine | stan-diagnose-helper.sh (bash, existing) |
| AI interpretation | Claude Haiku API via httr2 |
| Job status | Domino Jobs API + job-log.jsonl |
| Dependency graph | targets::tar_network() |
| Deployment | Domino App (single `run.sh` entry point) |

## Future Extensions

Not in MVP, but the design accommodates:

- **Job launching from UI** — requires auth confirmation flows, safety rails
- **Sclc project support** — add another Store; model hierarchy differs (dco_name dimension)
- **Multi-user** — would need auth layer; SQLite may need upgrading to PostgreSQL
- **Real-time fit monitoring** — tail CSV files during active sampling
- **Webhook notifications** — alert when severity changes from ok to bad
