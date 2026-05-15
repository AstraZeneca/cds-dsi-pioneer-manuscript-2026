# run_targets.sh
#
# General-purpose helper to run targets::tar_make with project-specific defaults.
# This is the base script used by sclc_targets.sh, pioneer_targets.sh, etc.
#
# Usage: ./run_targets.sh [OPTIONS] PROJECT_NAME
#
# Examples:
#   ./run_targets.sh sclc                 # run full sclc pipeline
#   ./run_targets.sh -c pioneer            # run pioneer with shortcut = TRUE
#   ./run_targets.sh -m '"a","b"' sclc    # run specific targets
#   ./run_targets.sh -d -m '"a"' sclc     # dry-run: print command only
#
# Function to display usage
usage() {
    echo "Usage: $0 [OPTIONS] PROJECT_NAME"
    echo ""
    echo "General-purpose targets runner. Use project-specific wrappers (sclc_targets.sh,"
    echo "pioneer_targets.sh) for convenience, or call this directly with a project name."
    echo ""
    echo "Options:"
    echo "  -s: Disable crew (use_crew = FALSE)"
    echo "  -c: Enable shortcut (shortcut = TRUE)"
    echo "  -d: Dry run (print command only)"
    echo "  -v: Never re-run targets (cue = tar_cue('never'))"
    echo "  -i: Invalidate specific targets"
    echo "  -m: Make specific targets (can be combined with -r)"
    echo "  -r: Target regex pattern(s) (can be combined with -m)"
    echo "  -n: Negate target regex (only affects -r patterns)"
    echo "  -D: Include downstream dependents (wrap selection in depends_on())"
    echo "  -b: Specify run/store name (default: main)"
    echo "  -p: Set SCLC_EXP_SUBDIR path for custom data directory"
    echo "  -u: Specify custom username (default: \$DOMINO_STARTING_USERNAME)"
    echo "  -l: Enable Laplace marginalization (sets ENABLE_LAPLACE=TRUE)"
    echo "  -k: Skip renv::restore()"
    echo ""
    echo "Arguments:"
    echo "  PROJECT_NAME: Name of targets project (e.g., sclc, pioneer)"
    exit 1
}

# Initialize variables
project_name=""
targets=""
use_crew="TRUE"
tar_run="main"
target_regex_negate="FALSE"
use_depends_on="FALSE"
shortcut="FALSE"
dry_run="FALSE"
use_cue_never="FALSE"
sclc_exp_subdir=""
skip_restore="FALSE"
custom_username=""
enable_laplace="FALSE"

# Parse command-line options
while getopts "i:m:r:b:p:u:h:sncdvklD" flag; do
    case "${flag}" in
        i) targets=${OPTARG};;
        m) make_targets=${OPTARG};;
        r) target_regex=${OPTARG};;
        n) target_regex_negate="TRUE";;
        D) use_depends_on="TRUE";;
        s) use_crew="FALSE";;
        c) shortcut="TRUE";;
        d) dry_run="TRUE";;
        v) use_cue_never="TRUE";;
        b) tar_run=${OPTARG};;
        p) sclc_exp_subdir=${OPTARG};;
        u) custom_username=${OPTARG};;
        l) enable_laplace="TRUE";;
        k) skip_restore="TRUE";;
        h) usage;;
        *) usage;;
    esac
done

# Shift the parsed options out of the argument list
shift $((OPTIND-1))

# Check if there's a remaining argument for project name
if [ $# -eq 1 ]; then
    project_name="$1"
elif [ $# -gt 1 ]; then
    echo "Error: Too many arguments, $#: $1" >&2
    usage
fi

# Project name is required
if [ -z "$project_name" ]; then
    echo "Error: Project name is required" >&2
    usage
fi

# Use custom username if provided, otherwise use DOMINO_STARTING_USERNAME
username="${custom_username:-$DOMINO_STARTING_USERNAME}"

target_store_dir="$DOMINO_DATASETS_DIR/analysis-results/$username/$project_name/$tar_run/_targets"

# Print diagnostic information about memory limits
echo "=========================================="
echo "Memory and Resource Diagnostics"
echo "=========================================="
echo "Hardware tier: ${DOMINO_HARDWARE_TIER_ID:-Not set}"
echo "Hostname: $(hostname)"
echo "Date: $(date)"
echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'N/A')"
echo "Git log: $(git log --oneline -1 2>/dev/null || echo 'N/A')"
echo ""

# Check cgroup memory limit (works for both cgroup v1 and v2)
echo "Cgroup Memory Limit:"
if [ -f /sys/fs/cgroup/memory/memory.limit_in_bytes ]; then
    CGROUP_LIMIT=$(cat /sys/fs/cgroup/memory/memory.limit_in_bytes)
    echo "  Raw: $CGROUP_LIMIT bytes"
    echo "  GB: $(python3 -c "print(f'{$CGROUP_LIMIT / 1024**3:.2f} GB')" 2>/dev/null || echo 'N/A')"
elif [ -f /sys/fs/cgroup/memory.max ]; then
    CGROUP_LIMIT=$(cat /sys/fs/cgroup/memory.max)
    echo "  Raw: $CGROUP_LIMIT"
    if [ "$CGROUP_LIMIT" != "max" ]; then
        echo "  GB: $(python3 -c "print(f'{$CGROUP_LIMIT / 1024**3:.2f} GB')" 2>/dev/null || echo 'N/A')"
    fi
else
    echo "  No cgroup memory limit found"
fi
echo ""

# Check system memory
echo "System Memory:"
free -h | head -n 2
echo ""

# Check ulimits
echo "Process Resource Limits (ulimit -a):"
ulimit -a | grep -E "memory|virtual|stack|data"
echo ""

echo "=========================================="
echo ""

# Run renv::restore() first (unless -k flag is set)
if [ "$skip_restore" = "FALSE" ]; then
    echo "Running renv::restore() to ensure packages are installed..."
    Rscript --vanilla -e "source('renv/activate.R'); renv::restore(prompt = FALSE)" 
    restore_status=$?

    if [ $restore_status -ne 0 ]; then
        echo "Error: renv::restore() failed with exit code $restore_status"
        exit $restore_status
    fi

    echo "renv::restore() completed successfully"
    echo ""
else
    echo "Skipping renv::restore() (-k flag set)"
    echo ""
fi

# Initialize the Rscript command
rscript_cmd="Rscript --no-restore --no-save" 

# Set SCLC_EXP_SUBDIR environment variable if provided
if [ -n "$sclc_exp_subdir" ]; then
    rscript_cmd+=" -e \"Sys.setenv(SCLC_EXP_SUBDIR = '$sclc_exp_subdir')\""
fi

# Enable Laplace marginalization if -l flag set
if [ "$enable_laplace" = "TRUE" ]; then
    rscript_cmd+=" -e \"Sys.setenv(ENABLE_LAPLACE = 'TRUE')\""
fi

# Add project name if provided
if [ -n "$project_name" ]; then
    rscript_cmd+=" -e \"Sys.setenv(TAR_PROJECT = '$project_name')\""
fi

# Set TAR_RUN so _targets.yaml store expression and tar_path_store() inside scripts reflect the correct store
rscript_cmd+=" -e \"Sys.setenv(TAR_RUN = '$tar_run')\""

# Process targets if provided
if [ -n "$targets" ]; then
    IFS=',' read -ra TARGET_ARRAY <<< "$targets"
    for target in "${TARGET_ARRAY[@]}"; do
        rscript_cmd+=" -e \"targets::tar_invalidate($target, store = '$target_store_dir')\""
    done
fi

# Set tar_cue if needed (must be before tar_make)
if [ "$use_cue_never" = "TRUE" ]; then
    rscript_cmd+=" -e \"targets::tar_option_set(cue = targets::tar_cue('never'))\""
fi

# Add the main tar_make command
TAR_MAKE_ARGS="callr_function = NULL, as_job = FALSE, use_crew = $use_crew, shortcut = $shortcut, reporter = 'terse', store = '$target_store_dir'"

# Build the targets expression combining both -m and -r options
targets_expr=""

# Add specific targets from -m option
if [ -n "$make_targets" ]; then
    targets_expr="$make_targets"
fi

# Add regex patterns from -r option
if [ -n "$target_regex" ]; then
    IFS=',' read -ra REGEX_ARRAY <<< "$target_regex"
    regex_expr=""
    for pat in "${REGEX_ARRAY[@]}"; do
        if [ -n "$regex_expr" ]; then
            regex_expr+=", matches('$pat', perl = TRUE)"
        else
            regex_expr+="matches('$pat', perl = TRUE)"
        fi
    done
    
    # Combine with make_targets if both are present
    if [ -n "$targets_expr" ]; then
        targets_expr+=", $regex_expr"
    else
        targets_expr="$regex_expr"
    fi
fi

# Execute tar_make with the appropriate targets
if [ -n "$targets_expr" ]; then
    # Build the inner selection expression
    if [ "$target_regex_negate" = "TRUE" ] && [ -n "$target_regex" ]; then
        # If negate is TRUE and we have regex patterns, negate only the regex part
        if [ -n "$make_targets" ]; then
            # Both -m and -r: combine specific targets with negated regex
            inner_expr="c($make_targets, !matches('$target_regex', perl = TRUE))"
        else
            # Only -r with negate
            inner_expr="!c($regex_expr)"
        fi
    else
        # Normal case: combine all targets
        inner_expr="c($targets_expr)"
    fi

    # Wrap in depends_on() if -D flag is set
    if [ "$use_depends_on" = "TRUE" ]; then
        final_expr="depends_on($inner_expr)"
    else
        final_expr="$inner_expr"
    fi

    rscript_cmd+=" -e \"targets::tar_make($final_expr, $TAR_MAKE_ARGS)\""
else
    # No specific targets or regex: run all
    rscript_cmd+=" -e \"targets::tar_make($TAR_MAKE_ARGS)\""
fi

echo "We are going to run: $rscript_cmd"

# Execute the Rscript command (or dry-run)
if [ "$dry_run" = "TRUE" ]; then
    # Just print and exit
    exit 0
else
    eval $rscript_cmd
fi