---
name: start-job
description: Start a Domino job with safety checks and intelligent target specification. Understands natural language target descriptions and constructs appropriate sclc_targets.sh commands.
user-invocable: true
allowed-tools: [Bash, AskUserQuestion, mcp__domino__start_job, mcp__domino__list_hardware_tiers, mcp__targets__tar_list_variations, mcp__targets__tar_manifest]
---

# Start Job Skill

Safely start a Domino job with pre-flight checks and interactive configuration.

## Purpose

This skill helps start Domino jobs while ensuring:
1. You're aware of any uncommitted git changes
2. You can specify targets using natural language
3. You **explicitly select a hardware tier (NO DEFAULTS - ALWAYS REQUIRED)**
4. You specify the correct pipeline store name
5. All configuration is confirmed before starting

## Key Concepts

**CRITICAL: Git Branch vs Pipeline Store Name**

These are two COMPLETELY SEPARATE things:

| Concept | What It Controls | How It's Specified | Example |
|---------|------------------|-------------------|---------|
| **Git Branch** | What CODE runs (model config, priors, Stan code) | Checked out in Domino workspace | `karim/process-noise`, `main` |
| **Pipeline Store Name** | WHERE results are saved | `-b` flag in `sclc_targets.sh` | `lfo`, `test`, `main`, `my-experiment` |

**Examples:**
- Code from branch `karim/process-noise` can save results to pipeline store `lfo`
- Code from branch `main` can save results to pipeline store `test`
- Pipeline store name is just an arbitrary label for organizing results

**Store Path:** `/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/{pipeline_store_name}/_targets`

## Usage

```
/start-job [pipeline_name] [target_spec]
```

**Arguments:**
- `pipeline_name` (optional) - Pipeline store name (e.g., "main", "test2", "karim/process-noise")
- `target_spec` (optional) - Natural language target specification

**Examples:**
```bash
/start-job test2 "lfo fits for ctdna aug"
/start-job main "posterior fits for no_oe apr"
/start-job "run tumor_ssls_lfo"              # Will ask for pipeline/model/dco
/start-job                                   # Full interactive mode
```

**Note:** Pipeline name corresponds to the store directory in `/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/{pipeline_name}/_targets`. This is passed to `sclc_targets.sh -b` flag. The `-b` flag specifies where targets data is stored, NOT which git branch to use. The git branch is determined by the Domino job's workspace (whatever branch is checked out in the Domino environment).

## Natural Language Target Specification

**Supported Keywords:**
- **Fit types**: "lfo", "prior", "posterior", "stan data", "stan_data"
- **Models**: Dynamically detected from sclc_targets.R (e.g., "ctdna", "no_oe", "no_ctdna", "ctdna_only", "all_trials", "no_covar", "no_hist", "no_pl")
- **DCOs**: Dynamically detected from sclc_targets.R (e.g., "apr", "aug")
- **Special**: "all" (run for all variations)

**Target Name Patterns:**
- Main fits: `tumor_ssls_{model}_{dco}` (e.g., `tumor_ssls_ctdna_aug`)
- Prior fits: `prior_tumor_ssls_{model}_{dco}` (e.g., `prior_tumor_ssls_no_oe_apr`)
- LFO fits: `tumor_ssls_lfo_{model}_{dco}` (e.g., `tumor_ssls_lfo_ctdna_aug`)
- Stan data: `tumor_ssls_stan_data_{type}_{model}_{dco}` (e.g., `tumor_ssls_stan_data_posterior_ctdna_aug`)

## Implementation Workflow

### 0. Parse Target Specification (if provided)

If user provides target specification in args or prompt:

**Step 0.1: Call tar_list_variations**
```bash
mcp__targets__tar_list_variations(project="/mnt/code")
```

Store returned variations and patterns for later use.

**Step 0.2: Extract Keywords from User Input**

Search for exact keyword matches (case-insensitive):
- Fit type: Check for "lfo", "prior", "posterior", "stan data"
- Model: Check for model names returned by tar_list_variations
- DCO: Check for dco names returned by tar_list_variations
- Quantifier: Check for "all" keyword

**Step 0.3: Determine Missing Variations**

Based on detected fit type:
- For "lfo": requires model + dco
- For "prior"/"posterior": requires model + dco
- For "stan data": requires fit_type (prior/posterior) + model + dco

**Step 0.4: Ask Questions for Missing Variations**

Use `AskUserQuestion` to fill gaps:

**Question: Pipeline Name** (if not provided as arg)
- **ASK DIRECTLY via text, not using AskUserQuestion tool**
- Simply ask: "What pipeline store name would you like to use?"
- User will respond with the name (e.g., "lfo", "test", "main", "my-experiment")
- Accept whatever name the user provides - this is just an arbitrary string for organizing results
- **IMPORTANT: Pipeline store name ≠ Git branch!**
  - **Pipeline store name** (`-b` flag): Where results are saved (e.g., `/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/lfo/_targets`)
  - **Git branch**: What code/model runs (determined by what's checked out in Domino workspace)
  - These are completely separate! Pipeline store "lfo" can run code from branch "karim/process-noise"

**Question: Model Selection** (if model not detected)
- Header: "Model"
- Question: "Which model variation should run?"
- Options (dynamically from tar_list_variations):
  - Models list (e.g., "ctdna", "no_oe", "no_ctdna", etc.)
  - "all - Run for all model variations"

**Question: DCO Selection** (if dco not detected)
- Header: "Data Cutoff"
- Question: "Which data cutoff should be used?"
- Options (dynamically from tar_list_variations):
  - DCO list (e.g., "apr", "aug")
  - "both - Run for both cutoffs"

**Question: Fit Type** (if needed for stan data)
- Header: "Fit Type"
- Question: "Prior or posterior fit?"
- Options:
  - "posterior - Fit to actual data"
  - "prior - Prior predictive sampling"

**Step 0.5: Construct Target Name(s)**

Use patterns from tar_list_variations:
- LFO: `tumor_ssls_lfo_{model}_{dco}`
- Posterior: `tumor_ssls_{model}_{dco}`
- Prior: `prior_tumor_ssls_{model}_{dco}`
- Stan data: `tumor_ssls_stan_data_{type}_{model}_{dco}`

If "all" selected for model, expand to comma-separated list:
```
tumor_ssls_lfo_ctdna_aug,tumor_ssls_lfo_no_oe_aug,tumor_ssls_lfo_no_ctdna_aug,...
```

**Step 0.6: Validate Target Exists**

Call `mcp__targets__tar_manifest(store="/mnt/data/analysis-results/$DOMINO_STARTING_USERNAME/sclc/main/_targets")`

Check if constructed target name(s) exist. If not found:
- Warn: "Target '{name}' not found in manifest"
- Ask: "Target doesn't exist in manifest. Continue anyway?" [Yes/No]
  - If No: return to model selection
  - If Yes: proceed (may be a new target being added)

### 1. Check Git Status

**Step 1.1: Check for Uncommitted Changes**

Use Bash to check for uncommitted changes:
```bash
git status --porcelain
```

**Step 1.2: Check for Unpushed Commits**

Check if local branch is ahead of remote:
```bash
git rev-list @{u}..HEAD --count 2>/dev/null || echo "0"
```

If count > 0, get the files changed in unpushed commits:
```bash
git diff --name-status @{u}..HEAD
```

**Step 1.3: Filter for Pipeline-Relevant Changes**

Only warn about changes that can affect the targets pipeline execution. **Ignore:**
- `.qmd` files (Quarto documentation)
- `plots/` directory and any plotting-related code
- `docs/` directory
- `.md` files (documentation)
- Changes to skill/agent definitions (`.claude/skills/`, `.claude/agents/`)

**Keep warnings for:**
- R code in `r/`, `targets/`
- Stan models in `stan/`
- Configuration files (`_targets.yaml`, `.Rprofile`, `renv.lock`)
- Shell scripts (`sclc_targets.sh`, etc.)

Example filtering:
```bash
# Get all uncommitted and unpushed files
all_changes=$(git status --porcelain; git diff --name-status @{u}..HEAD 2>/dev/null)

# Filter out non-pipeline files using grep
pipeline_changes=$(echo "$all_changes" | grep -v -E '\.(qmd|md)$|plots/|^docs/|\.claude/(skills|agents)/' | grep -E '\.(R|stan|sh|yaml|lock)$|^r/|^stan/|^targets/')
```

Only proceed with git warning if `pipeline_changes` is non-empty.

### 2. Get Current Branch

```bash
git branch --show-current
```

This shows what branch the user is currently on. Display this for informational purposes (helps user understand what CODE/model configuration will run in the job).

**CRITICAL: Git branch ≠ Pipeline store name**
- **Git branch**: Determines what CODE runs (model configuration, priors, Stan code, etc.)
- **Pipeline store name** (`-b` flag): Just a label for WHERE to save results
- Do NOT use git branch to auto-populate the pipeline store name - they are separate concepts!
- Example: Code from branch "karim/process-noise" can save results to pipeline store "lfo"

**Important:** The git branch shown is what's checked out locally. The Domino job will use whatever is checked out in its workspace, which should match if changes are pushed.

### 3. List Available Hardware Tiers

```bash
mcp__domino__list_hardware_tiers(project_id="67cf5d3f4abf2434ae3946a6")
```

Extract tier names, cores, memory, and GPUs for user selection.

### 4. Ask User Questions

Use `AskUserQuestion` to confirm configuration. Ask multiple questions:

**Question 1: Uncommitted/Unpushed Changes (if pipeline-relevant changes exist)**
- Header: "Git Status"
- Question: "You have uncommitted or unpushed changes that may affect the pipeline. Start job anyway?"
- Options:
  - "Yes, start with these changes"
  - "No, let me commit/push first"
- Show list of pipeline-relevant uncommitted and unpushed files in question description
- Note: Only show this question if there are pipeline-relevant changes (see Step 1.3 filtering)

**Question 2: Hardware Tier (REQUIRED - NO DEFAULTS)**
- Header: "Hardware"
- Question: "Which hardware tier should the job use?"
- **CRITICAL: Hardware tier MUST be explicitly specified - never use default tier**
- **IMPORTANT: Include ALL tiers returned from list_hardware_tiers as options**
- Format each option as: "tier_name - X cores, Y GB RAM" (add ", Z GPUs" if gpus > 0)
- Sort by type (CPU tiers first, then GPU tiers) and size
- Example option label: "ods-cpu-r6i (3XL) - 96 cores, 900 GB RAM"
- User can specify "Other" for custom tier ID
- **This question is mandatory and cannot be skipped**

**Question 3: Job Command**
- Header: "Command"
- Question: "What command should the job run?"
- Options:
  - If target(s) specified from step 0:
    - "sclc_targets.sh -m '{target_names}'" - Run specific target(s) (Recommended)
    - "sclc_targets.sh -r '{regex_pattern}'" - Use regex pattern instead
  - If no target specified:
    - "sclc_targets.sh" - Run full pipeline (Recommended)
  - User can specify "Other" for custom command

**Notes:**
- `-m "target1,target2"` runs specific comma-separated targets
- `-r "pattern"` runs targets matching regex (e.g., `"tumor_ssls_.*_aug"`)
- `-b pipeline_name` specifies pipeline store location (e.g., "main", "test2") - this is WHERE to store results, not which git branch to use
- `-s` disables crew for debugging
- The git branch used by the job is determined by what's checked out in the Domino workspace, not by the `-b` flag

### 5. Build Job Configuration

Based on user answers:
```
project_id: 67cf5d3f4abf2434ae3946a6 (SCLC project)
command: {base command}
hardware_tier_id: {selected tier ID - REQUIRED, must always be specified}
```

**CRITICAL:** `hardware_tier_id` is mandatory and must always be provided. Never omit this parameter or use a default.

**Note:** `environment_id` is typically not needed (uses project default).

### 6. Generate Job Title and Preview Command

**Generate Title:**
Agent should generate a descriptive title based on the job configuration:
- If specific target(s): "Run {target_name}" or "Run {N} targets: {short_description}"
- If full pipeline: "SCLC pipeline run - {pipeline_name}"
- Examples:
  - "Run tumor_ssls_lfo_ctdna_aug"
  - "Run LFO fits for ctdna (aug data)"
  - "Run 8 targets: all models LFO (aug)"
  - "SCLC pipeline run - main"

**Construct Full Command:**
Build complete sclc_targets.sh command with all flags:
```bash
sclc_targets.sh -b {pipeline_name} -m '{target_names}'
# or
sclc_targets.sh -b {pipeline_name} -r '{regex_pattern}'
# or
sclc_targets.sh -b {pipeline_name}  # full pipeline
```

**Show Preview:**
Display the full command that will be executed:
```
Job Configuration:
  Pipeline: {pipeline_name}
  Target(s): {target_names or "all targets"}
  Hardware: {tier_name} ({cores} cores, {memory} GB)
  Title: {generated_title}

Command that will run:
  sclc_targets.sh -b {pipeline_name} -m '{target_names}'
```

**Confirmation Question:**
- Header: "Confirm"
- Question: "Start job with this configuration?"
- Options:
  - "Yes, start job"
  - "No, cancel"

### 7. Start the Job

**CRITICAL:** Always provide `hardware_tier_id` - this parameter is mandatory and must never be omitted.

```bash
mcp__domino__start_job(
  project_id="67cf5d3f4abf2434ae3946a6",
  command="{full_command_with_flags}",
  title="{generated_title}",
  hardware_tier_id="{tier_id_from_user_selection}"  # REQUIRED - never omit
)
```

### 8. Report Results

Display:
- Job ID
- Job number
- Command
- Hardware tier
- Branch
- Link to job (if available in response)

**Success Format:**
```
✅ Job started successfully

Job ID: abc123xyz
Job Number: #456
Command: sclc_targets.sh
Branch: test2
Hardware: Large (8 cores, 32 GB RAM)
Title: SCLC run - test2

Monitor with: mcp__domino__get_job_status("abc123xyz")
```

**Error Format:**
```
❌ Failed to start job

Error: {error message}

Check:
- Project ID is correct
- Hardware tier ID is valid
- Command is valid
```

## Monitoring Job Progress

After starting a job, you can monitor its progress:
- Use the pipeline-monitor agent for interactive status updates
- Check Domino job status with `mcp__domino__get_job_status(job_id)`
- Check Domino job logs via the Domino UI

**Note:** The `_progress` file persists between runs and may show "dispatched" targets from previous runs. For accurate progress, check the Domino job logs directly.

## Safety Checks

### Uncommitted Changes Check
- Always run `git status --porcelain`
- If output is non-empty, warn user and ask for confirmation
- Show which files are modified/untracked
- User can cancel and commit first

### Branch Information
- Display current git branch for user reference
- The job will use whatever branch is checked out in the Domino workspace
- Ensure changes are pushed so Domino can access them

## Target Specification Examples

### Example 1: Complete specification
```
User: /start-job "run the lfo fits for ctdna aug"
Skill: Detects lfo + ctdna + aug
       Constructs: tumor_ssls_lfo_ctdna_aug
       Validates against manifest ✓
       Asks for hardware tier
       Command: sclc_targets.sh -b main -m 'tumor_ssls_lfo_ctdna_aug'
       Confirms and starts
```

### Example 2: Partial specification
```
User: /start-job "run tumor_ssls_lfo"
Skill: Detects lfo
       Missing: model, dco
       Asks: "Which model?" → User: "ctdna"
       Asks: "Which data cutoff?" → User: "aug"
       Constructs: tumor_ssls_lfo_ctdna_aug
       Validates and continues
```

### Example 3: Multiple targets
```
User: /start-job "run all lfo fits for aug"
Skill: Detects lfo + aug + "all"
       Expands to 8 targets (all models × aug)
       Asks: "Run 8 targets?"
       Command: sclc_targets.sh -b main -m 'tumor_ssls_lfo_ctdna_aug,tumor_ssls_lfo_no_oe_aug,...'
       Confirms and starts
```

### Example 4: Target not found
```
User: /start-job "run lfo for new_model aug"
Skill: Constructs: tumor_ssls_lfo_new_model_aug
       Validates: NOT in manifest ✗
       Warns: "Target not found"
       Asks: "Continue anyway?"
       User can cancel or proceed
```

## Common Scenarios

### Scenario 1: Quick start on current pipeline
```
User: /start-job
Assistant: Checks git status (clean)
          Shows current git branch: karim/process-noise
          (Note: This is what CODE will run, not where results are saved)
          Asks: "What pipeline store name would you like to use?"
          User responds: "main"
          (Note: Results will be saved to /mnt/data/analysis-results/.../sclc/main/_targets)
          Asks for hardware tier (shows all 13 available tiers)
          User selects: ods-cpu-m5 (L) - 20 cores, 160 GB RAM
          Suggests "sclc_targets.sh" command (full pipeline)
          Previews command: sclc_targets.sh -b main
          Auto-generates title: "SCLC pipeline run - main"
          User confirms
          Job starts with code from karim/process-noise branch, results saved to "main" store
```

### Scenario 2: Uncommitted changes
```
User: /start-job test2
Assistant: Detects 3 uncommitted files
          Shows list: M file1.R, M file2.stan, ?? file3.R
          Asks: "Start anyway?"
          User can cancel or proceed
```

### Scenario 3: Target-based job
```
User: /start-job main "lfo fits for ctdna aug"
Assistant: No uncommitted changes
          Pipeline: main
          Parses: lfo + ctdna + aug
          Constructs: tumor_ssls_lfo_ctdna_aug
          Asks for hardware tier (shows all 13 available tiers)
          User selects: ods-cpu-m5 (L) - 20 cores, 160 GB RAM
          Previews: sclc_targets.sh -b main -m 'tumor_ssls_lfo_ctdna_aug'
          Auto-generates title: "Run LFO fits for ctdna (aug data)"
          Starts job
```

## Project Configuration

**SCLC Project ID:** `67cf5d3f4abf2434ae3946a6`

This is hardcoded for the sclc project. If supporting multiple projects in the future, this should be parameterized.

## Error Handling

**Git not available:**
```
Warning: Could not check git status
Proceed without safety checks? [Yes/No]
```

**Hardware tiers API failure:**
```
Error: Could not list hardware tiers
Cannot proceed without hardware tier information.
Please check API connection or manually specify tier ID.
```

**Job start failure:**
```
Error: {API error message}
Common causes:
- Invalid project ID
- Invalid hardware tier ID
- Command syntax error
- Insufficient permissions
```

## Related Skills

- `/stan-flags` - Check model configuration flags before running
- Pipeline monitoring - Use pipeline-monitor agent after starting job

## Tips

1. **Before starting a job**, check what's already running:
   - Use pipeline-monitor agent: "what pipelines are running?"

2. **After starting a job**, monitor progress:
   - `mcp__domino__get_job_status("job_id")`
   - `mcp__domino__list_jobs(project_id, status="running")`

3. **Branch best practices**:
   - Commit changes before starting long-running jobs
   - Push branch to remote if needed for job access
   - Use descriptive job titles to identify runs later
