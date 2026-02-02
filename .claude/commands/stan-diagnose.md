# Stan/HMC Chain Diagnostics Skill

Use this skill to diagnose Stan MCMC sampling issues, particularly for complex hierarchical models.

## Usage

Invoke with `/stan-diagnose` followed by optional arguments:
- `/stan-diagnose` - Re-run last command (uses cached directory)
- `/stan-diagnose help` - Show available commands
- `/stan-diagnose <target> [on <store>]` - Quick health check of all chains
- `/stan-diagnose <pattern> [on <store>]` - Multi-target status (use wildcards like `*prior*ctdna*`)
- `/stan-diagnose <target> chain <N>` - Deep dive on specific chain (1-4)
- `/stan-diagnose <target> warmup` - Analyze warmup trajectories
- `/stan-diagnose <target> inits` - Compare initialization values across chains
- `/stan-diagnose <target> stuck` - Find iterations with extreme lp__
- `/stan-diagnose cache list` - Show cached diagnostics
- `/stan-diagnose cache clear` - Clear all cached data

Arguments:
- `<target>` - Name of the fit directory (e.g., `prior_tumor_ssls_ctdna_aug`)
- `<pattern>` - Glob pattern matching multiple targets (e.g., `*prior*ctdna_aug`, `tumor_ssls_*`)
- `on <store>` - Optional store name (e.g., `on main`, `on main3`). Defaults to `main`.

Flags:
- `-f`, `--force`, `--no-cache` - Bypass cache and force fresh read from CSV files

## Helper Script

The skill uses a bash helper script at `.claude/commands/stan-diagnose-helper.sh` for efficient file operations. The script handles:
- Dynamic header line counting (different Stan versions have different header sizes)
- Completed CSV files with footer comments
- Caching to avoid re-reading unchanged multi-GB files
- Large file handling via tail-based operations (no full file scans)

**To run diagnostics**, invoke the helper script via bash:
```bash
bash .claude/commands/stan-diagnose-helper.sh <output_dir> [command] [args]
```

## Multi-Target Pattern Support

When the user provides a pattern with wildcards (`*` or `?`), diagnose multiple targets:

1. Construct the glob path: `/mnt/data/analysis-results/karim_naguib/sclc/<store>/fit/<pattern>`
2. Find all matching directories
3. **Run status checks in PARALLEL** using multiple Bash calls in a single message
4. Present results with **run ID per target** (each target may have a different latest run)

**Parallel execution example** - When user asks for `*tumor*ctdna_aug on main`:
```
# In a SINGLE message, invoke multiple parallel bash calls:
bash .claude/commands/stan-diagnose-helper.sh "/mnt/data/.../prior_tumor_ssls_ctdna_aug" status
bash .claude/commands/stan-diagnose-helper.sh "/mnt/data/.../tumor_ssls_ctdna_aug" status
```

Example patterns:
- `*prior*ctdna*` - matches `prior_tumor_ssls_ctdna_aug`, `prior_tumor_ssls_no_ctdna_aug`, etc.
- `tumor_ssls_*_aug` - matches all aug variants
- `{prior_,}tumor_ssls_ctdna_aug` - matches exactly two targets (brace expansion)

## Caching System

The skill caches diagnostic output to avoid re-reading unchanged large CSV files (10-15GB each).

**Cache structure:**
```
.claude/.stan-diagnose-cache/
  <store>/
    <target>/
      <run_id>.cache
```

**Cache invalidation:** Uses file fingerprinting (mtime + size) for all 4 chain files. If any file changes, cache is invalidated.

**Cache commands:**
- `cache list` - Show all cached store/target/run combinations
- `cache clear` - Delete all cached data

**Force refresh:** Use `-f` flag to bypass cache: `/stan-diagnose prior_tumor_ssls_ctdna_aug -f`

## State Persistence

The last used output directory is cached in `.claude/.stan-diagnose-last-dir`. When `/stan-diagnose` is called with no arguments (or just a subcommand like `chain 2`), reuse the cached directory.

## Help Command

If the user runs `/stan-diagnose help`, display this summary:

```
Stan Chain Diagnostics - Available Commands:

  /stan-diagnose                           Re-run on last used directory
  /stan-diagnose help                      Show this help
  /stan-diagnose <target> [on <store>]     Quick status of all chains
  /stan-diagnose <pattern> [on <store>]    Multi-target status (parallel)
  /stan-diagnose <target> chain <N>        Deep dive on chain N (1-4)
  /stan-diagnose <target> warmup           Warmup progression analysis
  /stan-diagnose <target> inits            Compare initializations
  /stan-diagnose <target> stuck            Find iterations with lp__ < -100000
  /stan-diagnose cache list                Show cached diagnostics
  /stan-diagnose cache clear               Clear cache

Flags:
  -f, --force, --no-cache                  Bypass cache, force fresh read

Multi-target patterns (runs in parallel):
  /stan-diagnose *prior*ctdna* on main     All prior ctdna targets
  /stan-diagnose tumor_ssls_*_aug on main  All aug tumor_ssls variants

After first use, you can run subcommands without specifying the target:
  /stan-diagnose chain 2                   Deep dive on chain 2 (uses last dir)
  /stan-diagnose inits                     Compare inits (uses last dir)

Examples:
  /stan-diagnose prior_tumor_ssls_ctdna_aug on main
  /stan-diagnose tumor_ssls_aug chain 2
  /stan-diagnose -f                        (force refresh on last directory)
```

## Output Format

**Always present status as a markdown table with run ID per target:**

### target_name (Run: 202602022140)

| Chain | Progress | lp | stepsize | energy | divs (recent) | maxtree | Status |
|-------|----------|-----|----------|--------|---------------|---------|--------|
| 1 | 92% warmup | -82.1 | 0.032 | 1179 | 2 | 0 | ✓ healthy |
| 2 | 84% warmup | -56.2 | 0.025 | 1347 | 2 | 1 | ✓ healthy |
| 3 | 61% warmup | -66.7 | **0.007** | 1323 | **58** | 118 | ⚠️ adapting |
| 4 | 60% warmup | -122.2 | **0.004** | 1346 | **61** | 128 | ⚠️ adapting |

**Bold concerning values:**
- **stepsize < 0.01**: Small stepsize indicates difficult geometry
- **divs > 10**: High divergence count (especially during sampling)
- **maxtree > 50**: Frequent max treedepth hits in recent iterations
- **lp < -1e6**: Extreme negative log probability

## Implementation Notes

**Use the helper script** at `.claude/commands/stan-diagnose-helper.sh` for all file operations.

**Keep output formatting simple:**
- Collect raw numbers from script output
- Format in markdown table manually
- Bold concerning values based on thresholds above
- Present numbers, let human interpret

## Finding the Latest Run

The helper script automatically finds the latest run by:
1. Globbing `sf-ssm-log-space-2*.csv` (excludes profile files)
2. Extracting the 12-digit timestamp from filenames: `sf-ssm-log-space-YYYYMMDDHHMI-<chain>-<hash>.csv`

To manually find runs:
```bash
ls -t <output_dir>/sf-ssm-log-space-2*.csv | head -4
```

## Key Diagnostic Columns in Stan CSV Output

| Column | Name | Description | Healthy Range |
|--------|------|-------------|---------------|
| 1 | lp__ | Log posterior density | Higher = better; < -1e6 is bad |
| 2 | accept_stat__ | Acceptance probability | ~0.8 after adaptation |
| 3 | stepsize__ | Current step size | 0.01-0.1 typical; < 1e-10 is bad |
| 4 | treedepth__ | Tree depth | < 10 (max default) |
| 5 | n_leapfrog__ | Number of leapfrog steps | Varies |
| 6 | divergent__ | 1 if divergent transition | 0 is good; any 1s need investigation |
| 7 | energy__ | Hamiltonian energy | Similar across chains |
| 8+ | Model parameters | In declaration order | Model-dependent |

## Output Format

**Always present status as a markdown table with interpretation:**

| Chain | Progress | lp | stepsize | energy | divs (recent) | maxtree | Status |
|-------|----------|-----|----------|--------|---------------|---------|--------|
| 1 | 92% warmup | -82.1 | 0.032 | 1179 | 2 (+0) | 0 | ✓ healthy |
| 2 | 84% warmup | -56.2 | 0.025 | 1347 | 2 (+0) | 1 | ✓ healthy |
| 3 | 61% warmup | -66.7 | **0.007** | 1323 | **58** (+3) | 118 | ⚠️ adapting |
| 4 | 60% warmup | -122.2 | **0.004** | 1346 | **61** (+6) | 128 | ⚠️ adapting |

**Bold concerning values:**
- **stepsize < 0.01**: Small stepsize indicates difficult geometry
- **divs > 10**: High divergence count (especially during sampling)
- **recent divs > 0**: Active divergences in last 50 iterations
- **maxtree > 50**: Frequent max treedepth hits
- **lp < -1e6**: Extreme negative log probability

**Interpretation guidance:**
- Warmup divergences are less concerning than sampling divergences
- Chains recovering from bad inits will have high divs but should stabilize
- Energy should be similar across chains (within ~500)
- Stepsize should stabilize around 0.02-0.05 for most models

## Quick Health Check

**Use the helper script:**

```bash
bash .claude/commands/stan-diagnose-helper.sh <output_dir> status
```

The script outputs one line per chain with key diagnostics. Format this as a markdown table.

**Interpretation guidance:**
- Warmup divergences are less concerning than sampling divergences
- Chains recovering from bad inits will have high divs but should stabilize
- Energy should be similar across chains (within ~500)
- Stepsize should stabilize around 0.02-0.05 for most models

## Detecting Stuck Chains

Chains are stuck when:
- **lp__ < -1e6**: Extreme negative log probability indicates impossible parameter region
- **stepsize < 1e-10**: Stepsize collapsed trying to find valid moves
- **treedepth = 0 with divergent = 1**: Can't take any steps without diverging
- **No progress**: Line count not increasing over time

```bash
# Find iterations with extreme lp__
awk -F',' 'NR>48 && $1 < -100000 {print NR-48": lp="$1" step="$3" tree="$4" div="$6}' <file>.csv | head -20
```

## Analyzing Warmup Trajectories

When chains hit bad regions, examine:
1. **When they escaped**: First iteration with reasonable lp__
2. **What parameters changed**: Compare stuck vs escaped iterations
3. **Which parameters are problematic**: Usually hierarchical SDs or patient-level effects

```bash
# Find when chain escaped bad region (lp__ improved to > -1e6)
awk -F',' 'NR>48 && $1 > -1e6 && $1 < 0 {print NR-48": lp="$1; exit}' <file>.csv

# Show warmup progression every 50 iterations
awk -F',' 'NR>48 && (NR-48)%50==0 {printf "iter %d: lp=%.2f step=%.6f tree=%d div=%d energy=%.1f\n", NR-48, $1, $3, $4, $6, $7}' <file>.csv | head -15
```

## Multimodality Analysis

Multimodality occurs when the posterior has multiple distinct peaks (modes). This causes chains to converge to different regions, invalidating inference.

### Detecting Multimodality

**Compare sampling distributions across chains** (only valid after warmup):

```bash
# Compare lp__ and energy__ distributions across sampling chains
for i in 1 2 3 4; do
  FILE="sf-ssm-log-space-<RUN_ID>-$i-<HASH>.csv"
  adapted=$(grep -c "Adaptation terminated" "$FILE" 2>/dev/null || echo "0")
  if [ "$adapted" -eq 1 ]; then
    echo "Chain $i:"
    tail -n +549 "$FILE" | awk -F',' 'NR>1 {
      sum_lp+=$1; sum_lp2+=$1*$1;
      sum_e+=$7; sum_e2+=$7*$7;
      n++
    } END {
      if(n>0) {
        mean_lp=sum_lp/n; sd_lp=sqrt(sum_lp2/n - mean_lp*mean_lp);
        mean_e=sum_e/n; sd_e=sqrt(sum_e2/n - mean_e*mean_e);
        printf "  samples=%d, mean_lp=%.1f (sd=%.1f), mean_energy=%.1f (sd=%.1f)\n", n, mean_lp, sd_lp, mean_e, sd_e
      }
    }'
  else
    echo "Chain $i: still in warmup"
  fi
done
```

**Present results as a comparison table:**

| Chain | samples | mean lp (sd) | mean energy (sd) | Status |
|-------|---------|--------------|------------------|--------|
| 1 | 198 | -71.2 (56.0) | 1298.5 (174.6) | ✓ |
| 2 | 134 | -80.2 (43.4) | 1300.0 (203.4) | ✓ same mode |
| 3 | 89 | -75.1 (51.2) | 1305.2 (188.1) | ✓ same mode |
| 4 | 45 | **-250.3** (82.1) | **1580.2** (245.3) | ⚠️ **different mode?** |

**Signs of unimodality (good):**
- Mean lp values within ~1-2 standard deviations of each other
- Energy distributions nearly identical
- All chains exploring similar parameter ranges

**Signs of multimodality (concerning):**
- Mean lp differs by >2 SDs between chains
- Energy distributions significantly different (>500 difference in means)
- Chains show different parameter estimates for key quantities

### Model Design Concerns for Multimodality

**Common causes in hierarchical state-space models:**

| Concern | Mechanism | Risk Level | Mitigation |
|---------|-----------|------------|------------|
| **init/rate tradeoff** | Higher init + faster rate ≈ lower init + slower rate | Medium | Informative priors on init, sufficient early observations |
| **Process noise vs heterogeneity** | Both explain observation variance | Medium | Distinguish via temporal autocorrelation structure |
| **Label switching** | Mixture components can swap labels | High (if mixtures present) | Ordering constraints on mixture parameters |
| **Sign ambiguity** | Parameters can flip signs | Low (if rates on log scale) | Constrain signs via priors |
| **Weak identification** | Multiple parameter combos give similar likelihood | Medium | Stronger priors, more data |

**For tumor dynamics models specifically:**
- `log_decrease_rate = tr_loc + frac_log_decrease` - the tr/frac decomposition is identified because frac affects both decrease AND growth in opposite directions
- `states[t] = init + cumsum(rate * dt)` - init/rate tradeoff is a ridge, not multimodality, mitigated by hierarchical priors

### Distinguishing Bad Initialization vs True Multimodality

| Observation | Bad Init | True Multimodality |
|-------------|----------|-------------------|
| Chains differ during warmup | Expected | Expected |
| Chains converge after warmup | Yes - all to same region | No - stay in different regions |
| High warmup divergences | Common | Less common |
| lp gap persists post-warmup | No | Yes |
| Energy distributions match post-warmup | Yes | No |

**Key diagnostic**: If struggling chains converge to the same lp/energy region as healthy chains after warmup completes, it's initialization issues, not multimodality.

## Comparing Initializations

Check init JSON files (usually in /tmp/Rtmp*/):

```bash
# Find init files
ls /tmp/Rtmp*/init-*_*.json

# Compare key parameters across chains
for i in 1 2 3 4; do
  echo "=== Chain $i ==="
  python3 -c "
import json
d = json.load(open('<init_file>_$i.json'))
print('measure_sd:', d.get('measure_sd'))
print('tr_sd_patient_intercept:', d.get('tr_sd_patient_intercept'))
# Add other key parameters
"
done
```

## Common Problems and Solutions

### Problem: Chain stuck at lp = -1e50 with all-zero patient deviations
**Cause**: Fixed inits with zero heterogeneity - model can't explain patient-level data
**Solution**: Initialize patient raw effects with small spread: `rnorm(n_patients, sd=0.3)`

### Problem: Chains have different parameter counts than saved metric
**Cause**: Model flags changed (e.g., process noise enabled/disabled)
**Solution**: Disable metric_file or generate new metric from successful run with current config

### Problem: Chain 2 (or specific chain) always fails
**Cause**: Unlucky random initialization due to RNG state
**Solution**: Use fixed initializer with realistic patient heterogeneity

### Problem: Warmup extremely slow (hitting max treedepth)
**Cause**: Difficult posterior geometry, especially early in adaptation
**Solution**: Use pre-adapted metric from successful run, or be patient - this is sometimes normal

### Problem: High divergence count during warmup
**Cause**: Chains escaping bad initial region or adapting through difficult geometry
**Solution**: Monitor if divergences continue during sampling - warmup divergences are less critical

## Background Monitoring

The skill does NOT provide continuous monitoring. Each `/stan-diagnose` invocation is a one-shot status check.

**Options for monitoring:**

1. **Manual re-checks**: Ask for `/stan-diagnose` again when you want an update

2. **Parallel multi-target checks**: For multiple targets, multiple bash calls run in parallel and complete independently

3. **Background monitoring loop** (advanced): Start a bash loop that logs periodically:
```bash
while sleep 300; do
  echo "=== $(date '+%H:%M:%S') ===" >> /tmp/stan-monitor.log
  bash .claude/commands/stan-diagnose-helper.sh <output_dir> status >> /tmp/stan-monitor.log
done
```

The caching system ensures repeated checks are fast when files haven't changed.

## Efficient File Reading

Stan CSV files can be huge (10-15GB with ~850 lines, each line ~15MB). The helper script handles this efficiently:

- **Dynamic header counting**: `grep -c "^#"` to handle different Stan versions
- **Tail-based last line**: `tail -10 | grep -v "^#" | tail -1` to skip footer comments
- **Recent divergence counting**: `tail -55 | grep -v "^#" | tail -50` for last 50 data rows

**Never use** `cat`, `head -n <large>`, or full-file `grep` on these files.

For manual queries on completed files:
```bash
# Last data line (skip footer comments)
tail -10 file.csv | grep -v "^#" | tail -1 | cut -d',' -f1,3,6,7

# Use duckdb for complex queries
duckdb -c "SELECT \"lp__\", \"stepsize__\", \"divergent__\" FROM read_csv('file.csv', skip=51, header=true, max_line_size=100000000) LIMIT 10"
```
