# Stan/HMC Chain Diagnostics Skill

Use this skill to diagnose Stan MCMC sampling issues, particularly for complex hierarchical models.

## Usage

Invoke with `/stan-diagnose` followed by optional arguments:
- `/stan-diagnose` or `/stan-diagnose status` - Quick health check of all chains
- `/stan-diagnose chain <N>` - Deep dive on specific chain
- `/stan-diagnose warmup` - Analyze warmup trajectories for problematic geometry
- `/stan-diagnose inits` - Compare initialization values across chains
- `/stan-diagnose stuck` - Find iterations where chains hit bad regions

## Implementation Notes

**Keep it simple:**
- Use basic bash arithmetic `$((x * 100 / y))` - NO bc required
- Collect raw numbers first, format in markdown table manually
- Don't try to automate bold text or status emoji in bash
- Use `grep -c`, `wc -l`, `cut`, `tail` - all standard tools
- Present numbers, let human interpret against thresholds

## Finding the Latest Run

First, identify the output directory and latest run:

```bash
# Find most recent CSV files
ls -lt <output_dir>/sf-ssm-log-space-*.csv | head -8

# Extract run timestamp from filename pattern: sf-ssm-log-space-YYYYMMDDHHMI-<chain>-<hash>.csv
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

**Simple diagnostic script (no dependencies):**

```bash
# Navigate to fit directory
cd <output_dir>

# Quick status for all chains
for i in 1 2 3 4; do
  f=$(ls sf-ssm-log-space-<RUN_ID>-${i}-*.csv 2>/dev/null | head -1)
  if [ -f "$f" ]; then
    lines=$(wc -l < "$f")
    data_lines=$((lines - 49))

    # Get last sample values
    last=$(tail -1 "$f")
    lp=$(echo "$last" | cut -d',' -f1)
    step=$(echo "$last" | cut -d',' -f3)
    energy=$(echo "$last" | cut -d',' -f7)

    # Count diagnostics
    divs=$(tail -n +49 "$f" | cut -d',' -f6 | grep -c "^1" 2>/dev/null || echo "0")
    divs_recent=$(tail -50 "$f" | cut -d',' -f6 | grep -c "^1" 2>/dev/null || echo "0")
    maxtree=$(tail -n +49 "$f" | cut -d',' -f4 | grep -c "^10" 2>/dev/null || echo "0")

    # Check if adapted
    if grep -q "Adaptation terminated" "$f"; then
      pct=$((data_lines * 100 / 1000))
      stage="${pct}% sampling"
    else
      pct=$((data_lines * 100 / 500))
      stage="${pct}% warmup"
    fi

    printf "Chain %d: %s | lp=%.1f | step=%.4f | energy=%.0f | divs=%d (+%d) | maxtree=%d\n" \
      $i "$stage" "$lp" "$step" "$energy" $divs $divs_recent $maxtree
  fi
done
```

**Present as markdown table:**

| Chain | Progress | lp | stepsize | energy | divs (recent) | maxtree | Status |
|-------|----------|-----|----------|--------|---------------|---------|--------|
| 1 | 97% warmup | 115.1 | 0.021 | 1127 | 68 (+9) | 131 | ⚠️ adapting |
| 2 | 99% warmup | 125.1 | 0.0089 | 1083 | 66 (+9) | 130 | ⚠️ adapting |
| 3 | 93% warmup | 53.5 | 0.029 | 1299 | 69 (+6) | 142 | ⚠️ adapting |
| 4 | 95% warmup | -142.8 | 0.0077 | 1436 | 68 (+8) | 136 | ⚠️ adapting |

**Manual interpretation (don't automate):**
- Compare values to thresholds in "Bold concerning values" section above
- Assess energy convergence (all chains within ~500 of each other = good)
- Check if warmup divergences stabilize during sampling
- Note: Avoid using `bc` or complex conditionals - just present numbers

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

**Option 1: Simple periodic check**

```bash
# Check every 60 seconds
while sleep 60; do
  echo "=== $(date '+%H:%M:%S') ==="
  cd <output_dir>
  for i in 1 2 3 4; do
    f=$(ls sf-ssm-log-space-<RUN_ID>-${i}-*.csv 2>/dev/null | head -1)
    [ ! -f "$f" ] && continue

    lines=$(wc -l < "$f")
    last=$(tail -1 "$f")
    lp=$(echo "$last" | cut -d',' -f1)
    step=$(echo "$last" | cut -d',' -f3)

    # Simple progress
    data=$((lines - 49))
    pct=$((data * 100 / 500))

    printf "Chain %d: %3d%% | lp=%8.1f | step=%.4f\n" $i $pct "$lp" "$step"
  done
  echo ""
done
```

**Option 2: Run diagnostics script in background**

Use `run_in_background=true` for bash commands to monitor while working on other tasks.

## Efficient File Reading

Stan CSV files can be huge (multi-GB). Avoid `cat`, `head -n`, or loading entire file:

```bash
# Use awk for filtered extraction
awk -F',' 'NR>48 && NR<100 {print $1,$3,$4,$6,$7}' file.csv

# Use cut for specific columns from last line
tail -1 file.csv | cut -d',' -f1,3,4,6,7

# Use duckdb for complex queries (set max_line_size for wide files)
duckdb -c "SELECT \"lp__\", \"stepsize__\", \"divergent__\", \"energy__\" FROM read_csv('file.csv', skip=47, header=true, max_line_size=100000000) LIMIT 10"
```
