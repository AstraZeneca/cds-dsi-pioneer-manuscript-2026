#!/bin/bash
# Stan Chain Diagnostics Helper Script
# Usage: ./stan-diagnose-helper.sh <output_dir> [command] [args]
#
# Commands:
#   status       - Quick health check of all chains (default)
#   chain <N>    - Deep dive on specific chain (1-4)
#   warmup       - Show warmup progression for all chains
#   stuck        - Find iterations with extreme lp__
#   inits        - Compare initialization values across chains
#   monitor      - Start background monitoring (logs every 5 min)
#   help         - Show this help message

set -e

# Handle flags anywhere in args
FORCE_REFRESH=false
for arg in "$@"; do
    if [[ "$arg" == "-f" || "$arg" == "--force" || "$arg" == "--no-cache" ]]; then
        FORCE_REFRESH=true
    fi
    if [[ "$arg" == "-h" || "$arg" == "--help" || "$arg" == "help" ]]; then
        echo "Stan Chain Diagnostics Helper Script"
        echo ""
        echo "Usage: $0 <output_dir> [command] [args]"
        echo "       $0 [command] [args]              (uses cached directory)"
        echo ""
        echo "Commands:"
        echo "  status         Quick health check of all chains (default)"
        echo "  chain <N>      Deep dive on specific chain (1-4)"
        echo "  warmup         Show warmup progression for all chains"
        echo "  stuck          Find iterations with lp__ < -100000"
        echo "  inits          Compare initialization values across chains"
        echo "  monitor        Start background monitoring (logs every 5 min)"
        echo "  cache [list|clear]  Manage output cache"
        echo "  help           Show this help message"
        echo ""
        echo "Flags:"
        echo "  -f, --force, --no-cache   Bypass cache and force fresh read"
        echo ""
        echo "Examples:"
        echo "  $0 /path/to/fit/dir                    # Quick status check"
        echo "  $0 /path/to/fit/dir chain 2            # Deep dive on chain 2"
        echo "  $0 chain 2                             # Uses last directory"
        echo "  $0                                     # Re-run status on last dir"
        echo ""
        echo "Cache files:"
        echo "  Directory: .claude/.stan-diagnose-last-dir"
        echo "  Output:    .claude/.stan-diagnose-output-cache (mtime+size based)"
        echo ""
        echo "Stan CSV columns: 1=lp__, 2=accept_stat__, 3=stepsize__,"
        echo "                  4=treedepth__, 5=n_leapfrog__, 6=divergent__, 7=energy__"
        exit 0
    fi
done

# Use project-local cache if .claude dir exists, otherwise /tmp
if [ -d ".claude" ]; then
    CACHE_DIR=".claude/.stan-diagnose-cache"
    LAST_DIR_FILE=".claude/.stan-diagnose-last-dir"
else
    CACHE_DIR="/tmp/.stan-diagnose-cache"
    LAST_DIR_FILE="/tmp/.stan-diagnose-last-dir"
fi

# Extract store and target from output directory path
# Expected format: .../sclc/<store>/fit/<target>/
extract_store_target() {
    local dir="$1"
    # Extract store (parent of fit directory)
    STORE=$(echo "$dir" | grep -oP 'sclc/\K[^/]+(?=/fit)' || echo "unknown")
    # Extract target (directory name after fit/)
    TARGET=$(basename "$dir")
}

# Get file fingerprint (mtime + size) - fast way to detect changes
get_fingerprint() {
    local file="$1"
    if [ -f "$file" ]; then
        stat -c "%Y-%s" "$file" 2>/dev/null || echo "0-0"
    else
        echo "0-0"
    fi
}

# Get fingerprints for all chain files in a run
get_run_fingerprint() {
    local dir="$1"
    local run_id="$2"
    local fp=""
    for i in 1 2 3 4; do
        local f=$(ls "$dir"/sf-ssm-log-space-${run_id}-${i}-*.csv 2>/dev/null | head -1)
        fp="${fp}:$(get_fingerprint "$f")"
    done
    echo "$fp"
}

# Get cache file path for a specific store/target/run
get_cache_path() {
    local store="$1"
    local target="$2"
    local run_id="$3"
    echo "${CACHE_DIR}/${store}/${target}/${run_id}.cache"
}

# Check if cached output is still valid
check_cache() {
    local dir="$1"
    local run_id="$2"

    # Force refresh bypasses cache
    if [ "$FORCE_REFRESH" = true ]; then
        return 1
    fi

    extract_store_target "$dir"
    local cache_file=$(get_cache_path "$STORE" "$TARGET" "$run_id")

    if [ ! -f "$cache_file" ]; then
        return 1  # No cache
    fi

    # Read cached fingerprint (first line)
    local cached_fp=$(head -1 "$cache_file")

    # Check if files changed
    local current_fp=$(get_run_fingerprint "$dir" "$run_id")
    if [ "$cached_fp" != "$current_fp" ]; then
        return 1
    fi

    return 0  # Cache valid
}

# Save output to cache
save_cache() {
    local dir="$1"
    local run_id="$2"
    local output="$3"

    extract_store_target "$dir"
    local cache_file=$(get_cache_path "$STORE" "$TARGET" "$run_id")
    local cache_dir=$(dirname "$cache_file")

    # Create cache directory structure
    mkdir -p "$cache_dir"

    local fp=$(get_run_fingerprint "$dir" "$run_id")
    {
        echo "$fp"
        echo "$output"
    } > "$cache_file"
}

# Get cached output (lines 2+)
get_cached_output() {
    local dir="$1"
    local run_id="$2"

    extract_store_target "$dir"
    local cache_file=$(get_cache_path "$STORE" "$TARGET" "$run_id")
    tail -n +2 "$cache_file"
}

# List all cached runs
list_caches() {
    echo "Cached diagnostics:"
    if [ -d "$CACHE_DIR" ]; then
        find "$CACHE_DIR" -name "*.cache" -type f | while read -r f; do
            rel_path="${f#$CACHE_DIR/}"
            store=$(echo "$rel_path" | cut -d'/' -f1)
            target=$(echo "$rel_path" | cut -d'/' -f2)
            run_id=$(basename "$f" .cache)
            fp=$(head -1 "$f")
            echo "  $store/$target/$run_id"
        done
    else
        echo "  (no caches)"
    fi
}

# Known subcommands
SUBCOMMANDS="status|chain|warmup|stuck|inits|monitor|cache"

# Filter out flags from positional arguments
POSITIONAL_ARGS=()
for arg in "$@"; do
    if [[ "$arg" != "-f" && "$arg" != "--force" && "$arg" != "--no-cache" ]]; then
        POSITIONAL_ARGS+=("$arg")
    fi
done

ARG1="${POSITIONAL_ARGS[0]:-}"
ARG2="${POSITIONAL_ARGS[1]:-}"
ARG3="${POSITIONAL_ARGS[2]:-}"

# Determine if first arg is a directory or subcommand
USE_CACHE=false
if [ -z "$ARG1" ] || [ "$ARG1" = "." ]; then
    USE_CACHE=true
    COMMAND="${ARG2:-status}"
    CHAIN_NUM="${ARG3:-}"
elif [[ "$ARG1" =~ ^($SUBCOMMANDS)$ ]]; then
    USE_CACHE=true
    COMMAND="$ARG1"
    CHAIN_NUM="${ARG2:-}"
else
    OUTPUT_DIR="$ARG1"
    COMMAND="${ARG2:-status}"
    CHAIN_NUM="${ARG3:-}"
fi

if [ "$USE_CACHE" = true ]; then
    if [ -f "$LAST_DIR_FILE" ]; then
        OUTPUT_DIR=$(cat "$LAST_DIR_FILE")
    else
        echo "No cached directory found. Please specify a directory."
        echo "Usage: $0 <output_dir> [command] [args]"
        exit 1
    fi
fi

# Check if OUTPUT_DIR contains wildcards (multi-target pattern)
MULTI_TARGET=false
if [[ "$OUTPUT_DIR" == *"*"* ]] || [[ "$OUTPUT_DIR" == *"?"* ]]; then
    MULTI_TARGET=true
fi

# Count header lines dynamically (comment lines starting with # plus column header)
# PERFORMANCE: Only reads first 200KB - headers are always at the top of the file
count_header_lines() {
    local file="$1"
    local comment_lines=$(head -c 200000 "$file" 2>/dev/null | grep -c "^#" || echo "0")
    echo $((comment_lines + 1))  # +1 for column names row
}

# Estimate total data lines without full-file scan
# Uses file size and sample line length from the last line
estimate_line_count() {
    local file="$1"
    local header_lines="$2"
    local file_size=$(stat -c%s "$file" 2>/dev/null || echo "0")
    # Get a sample line length from the last data line (tail is fast - seeks from end)
    local sample_line=$(tail -c 60000000 "$file" | grep -v "^#" | tail -1)
    local line_len=${#sample_line}
    if [ "$line_len" -gt 0 ]; then
        # +1 for newline character
        local total_lines=$(( file_size / (line_len + 1) ))
        local data_lines=$(( total_lines - header_lines ))
        echo "$data_lines"
    else
        echo "0"
    fi
}

# Check adaptation status: infers from estimated line count vs num_warmup
# With ~15MB lines, the "Adaptation terminated" comment is buried deep in multi-GB files
# and can't be found via tail. Instead, infer from data: if lines > warmup, we're sampling.
# Args: file, header_lines, num_warmup
# Returns: 0 if adapted (sampling), 1 if still in warmup
check_adapted() {
    local file="$1"
    local header_lines="$2"
    local num_warmup="$3"
    local data_lines=$(estimate_line_count "$file" "$header_lines")
    # Add 10% margin to account for estimation error
    local threshold=$(( num_warmup + num_warmup / 10 ))
    [ "$data_lines" -gt "$threshold" ]
}

# Read config values from Stan CSV header (first ~200KB)
# These are comment lines like "# num_warmup = 300"
read_header_config() {
    local file="$1"
    local key="$2"
    head -c 200000 "$file" 2>/dev/null | grep "$key" | head -1 | grep -oP '\d+' || echo ""
}

# Find the latest run (exclude profile files)
find_latest_run() {
    local dir="$1"
    ls -t "$dir"/sf-ssm-log-space-2*.csv 2>/dev/null | head -1 | grep -oP '\d{12}' | head -1
}

# Function to run status on a single directory
run_single_status() {
    local dir="$1"
    local run_id=$(find_latest_run "$dir")

    if [ -z "$run_id" ]; then
        echo "  (no Stan output files)"
        return
    fi

    extract_store_target "$dir"

    # Check cache
    if [ "$FORCE_REFRESH" != true ] && check_cache "$dir" "$run_id"; then
        echo "  (cached)"
        get_cached_output "$dir" "$run_id" | sed 's/^/  /'
        return
    fi

    local output=""
    for i in 1 2 3 4; do
        local FILE=$(ls "$dir"/sf-ssm-log-space-${run_id}-${i}-*.csv 2>/dev/null | head -1)
        if [ -f "$FILE" ]; then
            local header_lines=$(count_header_lines "$FILE")
            local num_warmup=$(read_header_config "$FILE" "num_warmup")
            num_warmup=${num_warmup:-300}
            local num_samples=$(read_header_config "$FILE" "num_samples")
            num_samples=${num_samples:-500}
            local last_line=$(tail -c 60000000 "$FILE" | grep -v "^#" | tail -1)
            local lp=$(echo "$last_line" | cut -d',' -f1)
            local stepsize=$(echo "$last_line" | cut -d',' -f3)
            local energy=$(echo "$last_line" | cut -d',' -f7)
            local adapted=0
            if check_adapted "$FILE" "$header_lines" "$num_warmup"; then adapted=1; fi
            local divs_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f6 | grep -c "^1" 2>/dev/null | tr -d '\n' || echo "0")

            local data_lines=$(estimate_line_count "$FILE" "$header_lines")
            local stage pct
            if [ "$adapted" -eq 1 ]; then
                local sampling_lines=$((data_lines - num_warmup))
                [ "$sampling_lines" -lt 0 ] && sampling_lines=0
                pct=$((sampling_lines * 100 / num_samples))
                stage="samp"
                data_lines=$sampling_lines
            else
                pct=$((data_lines * 100 / num_warmup))
                stage="warm"
            fi

            local line_out="  Ch$i: ${pct}% ${stage}, lp=${lp}, div=${divs_recent}"
            echo "$line_out"
            output="${output}${line_out}"$'\n'
        fi
    done

    # Save to cache
    save_cache "$dir" "$run_id" "$output"
}

# Handle multi-target pattern
if [ "$MULTI_TARGET" = true ]; then
    echo "=== Multi-target Status ==="
    echo "Pattern: $OUTPUT_DIR"
    echo ""

    # Expand the glob pattern
    matched_dirs=()
    for d in $OUTPUT_DIR; do
        if [ -d "$d" ]; then
            matched_dirs+=("$d")
        fi
    done

    if [ ${#matched_dirs[@]} -eq 0 ]; then
        echo "No matching directories found"
        exit 1
    fi

    echo "Found ${#matched_dirs[@]} matching target(s):"
    echo ""

    for dir in "${matched_dirs[@]}"; do
        target=$(basename "$dir")
        run_id=$(find_latest_run "$dir")
        echo "--- $target (run: ${run_id:-none}) ---"
        if [ -n "$run_id" ]; then
            run_single_status "$dir"
        else
            echo "  (no Stan output)"
        fi
        echo ""
    done

    exit 0
fi

# Single target mode (original behavior)
RUN_ID=$(find_latest_run "$OUTPUT_DIR")
if [ -z "$RUN_ID" ]; then
    echo "No Stan output files found in $OUTPUT_DIR"
    exit 1
fi

# Cache the directory for next time
echo "$OUTPUT_DIR" > "$LAST_DIR_FILE"

echo "Run ID: $RUN_ID"
echo "Directory: $OUTPUT_DIR"
echo "---"

# Stan CSV columns:
# 1=lp__, 2=accept_stat__, 3=stepsize__, 4=treedepth__, 5=n_leapfrog__, 6=divergent__, 7=energy__

case "$COMMAND" in
    status)
        # Check if we can use cached output
        if check_cache "$OUTPUT_DIR" "$RUN_ID"; then
            echo "(cached - use -f to force refresh)"
            get_cached_output "$OUTPUT_DIR" "$RUN_ID"
        else
            # Compute fresh output
            output=""
            for i in 1 2 3 4; do
                FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
                if [ -f "$FILE" ]; then
                    header_lines=$(count_header_lines "$FILE")

                    # Get config from header (reads only first 200KB)
                    num_warmup=$(read_header_config "$FILE" "num_warmup")
                    num_warmup=${num_warmup:-300}
                    num_samples=$(read_header_config "$FILE" "num_samples")
                    num_samples=${num_samples:-500}

                    # Get last non-comment line values efficiently (tail seeks from end)
                    last_line=$(tail -c 60000000 "$FILE" | grep -v "^#" | tail -1)
                    lp=$(echo "$last_line" | cut -d',' -f1)
                    stepsize=$(echo "$last_line" | cut -d',' -f3)
                    treedepth=$(echo "$last_line" | cut -d',' -f4)
                    energy=$(echo "$last_line" | cut -d',' -f7)

                    # Check adaptation status (inferred from line count vs warmup)
                    adapted=0
                    if check_adapted "$FILE" "$header_lines" "$num_warmup"; then adapted=1; fi

                    # Recent divergences (last 50 data lines, exclude comments)
                    divs_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f6 | grep -c "^1" 2>/dev/null | tr -d '\n' || echo "0")
                    maxtree_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f4 | grep -c "^10" 2>/dev/null | tr -d '\n' || echo "0")

                    # Estimate data lines from file size (avoids wc -l full scan)
                    data_lines=$(estimate_line_count "$FILE" "$header_lines")

                    # Calculate progress
                    if [ "$adapted" -eq 1 ]; then
                        sampling_lines=$((data_lines - num_warmup))
                        [ "$sampling_lines" -lt 0 ] && sampling_lines=0
                        pct=$((sampling_lines * 100 / num_samples))
                        stage="sampling"
                        data_lines=$sampling_lines
                    else
                        pct=$((data_lines * 100 / num_warmup))
                        stage="warmup"
                    fi

                    # Warn if problematic
                    warn=""
                    if [[ "$lp" == *"e+"* ]] && [[ "$lp" == "-"* ]]; then
                        warn=" [EXTREME LP!]"
                    fi
                    if [[ "$stepsize" == *"e-"* ]]; then
                        exp=$(echo "$stepsize" | grep -oP 'e-\d+' | grep -oP '\d+')
                        if [ -n "$exp" ] && [ "$exp" -gt 10 ]; then
                            warn="$warn [TINY STEPSIZE!]"
                        fi
                    fi

                    line_out="Chain $i: ~${data_lines} ${stage} (~${pct}%), lp=$lp, step=$stepsize, energy=$energy, divs_recent=$divs_recent, maxtree_recent=$maxtree_recent, adapted=$adapted$warn"
                else
                    line_out="Chain $i: FILE NOT FOUND"
                fi
                echo "$line_out"
                output="${output}${line_out}"$'\n'
            done
            # Save to cache
            save_cache "$OUTPUT_DIR" "$RUN_ID" "$output"
        fi
        ;;

    stuck)
        echo "Finding iterations with lp__ < -100000..."
        echo "Note: Only checks last 200 iterations (use chain deep dive for full history)"
        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                echo "=== Chain $i ==="
                # Only check recent iterations to avoid full-file scan
                tail -c 5500000000 "$FILE" | grep -v "^#" | tail -200 | awk -F',' '$1 < -100000 {print NR": lp="$1" step="$3" tree="$4" div="$6" energy="$7}' | head -10
            fi
        done
        ;;

    warmup)
        echo "Warmup progression (every 50 iterations)..."
        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                header_lines=$(count_header_lines "$FILE")
                echo "=== Chain $i ==="
                # Use head to limit input to ~800 lines (enough for warmup), then awk
                head -800 "$FILE" | awk -F',' -v hdr="$header_lines" 'NR>hdr && (NR-hdr)%50==0 {printf "iter %d: lp=%.2f step=%.6f tree=%d div=%d energy=%.1f\n", NR-hdr, $1, $3, $4, $6, $7}' | head -15
            fi
        done
        ;;

    monitor)
        LOG="/tmp/chain_monitor_${RUN_ID}.log"
        echo "Starting background monitor, logging to $LOG"
        echo "Press Ctrl+C to stop"
        while true; do
            echo "=== $(date) ===" >> "$LOG"
            for i in 1 2 3 4; do
                FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
                if [ -f "$FILE" ]; then
                    header_lines=$(count_header_lines "$FILE")
                    num_warmup=$(read_header_config "$FILE" "num_warmup")
                    num_warmup=${num_warmup:-300}
                    num_samples=$(read_header_config "$FILE" "num_samples")
                    num_samples=${num_samples:-500}
                    last_line=$(tail -c 60000000 "$FILE" | grep -v "^#" | tail -1)
                    lp=$(echo "$last_line" | cut -d',' -f1)
                    stepsize=$(echo "$last_line" | cut -d',' -f3)
                    energy=$(echo "$last_line" | cut -d',' -f7)
                    divs_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f6 | grep -c "^1" 2>/dev/null || echo "0")
                    maxtree_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f4 | grep -c "^10" 2>/dev/null || echo "0")
                    adapted=0
                    if check_adapted "$FILE" "$header_lines" "$num_warmup"; then adapted=1; fi
                    data_lines=$(estimate_line_count "$FILE" "$header_lines")
                    if [ "$adapted" -eq 1 ]; then
                        sampling_lines=$((data_lines - num_warmup))
                        [ "$sampling_lines" -lt 0 ] && sampling_lines=0
                        pct=$((sampling_lines * 100 / num_samples))
                        stage="~$sampling_lines sampling"
                    else
                        pct=$((data_lines * 100 / num_warmup))
                        stage="~$data_lines warmup"
                    fi
                    echo "Chain $i: (~$data_lines data, $stage, ~$pct%) lp=$lp step=$stepsize energy=$energy divs_recent=$divs_recent maxtree_recent=$maxtree_recent adapted=$adapted" >> "$LOG"
                fi
            done
            echo "Logged at $(date)"
            sleep 300
        done
        ;;

    chain)
        if [ -z "$CHAIN_NUM" ] || ! [[ "$CHAIN_NUM" =~ ^[1-4]$ ]]; then
            echo "Error: Please specify chain number (1-4)"
            echo "Usage: $0 <output_dir> chain <N>"
            exit 1
        fi
        FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${CHAIN_NUM}-*.csv 2>/dev/null | head -1)
        if [ ! -f "$FILE" ]; then
            echo "Chain $CHAIN_NUM: FILE NOT FOUND"
            exit 1
        fi

        echo "=== Chain $CHAIN_NUM Deep Dive ==="
        echo "File: $FILE"
        echo ""

        header_lines=$(count_header_lines "$FILE")
        num_warmup=$(read_header_config "$FILE" "num_warmup")
        num_warmup=${num_warmup:-300}
        num_samples=$(read_header_config "$FILE" "num_samples")
        num_samples=${num_samples:-500}
        adapted=0
        if check_adapted "$FILE" "$header_lines" "$num_warmup"; then adapted=1; fi
        seed=$(head -c 200000 "$FILE" | grep "^# *seed" | head -1 | grep -oP '\d+' || echo "unknown")
        init_file=$(head -c 200000 "$FILE" | grep "^# init = " | head -1 | sed 's/.*= //')

        echo "Configuration:"
        echo "  Seed: $seed"
        echo "  Warmup: $num_warmup, Samples: $num_samples"
        echo "  Init file: $init_file"
        echo ""

        # Progress (estimated from file size)
        data_lines=$(estimate_line_count "$FILE" "$header_lines")
        if [ "$adapted" -eq 1 ]; then
            sampling_lines=$((data_lines - num_warmup))
            [ "$sampling_lines" -lt 0 ] && sampling_lines=0
            pct=$((sampling_lines * 100 / num_samples))
            echo "Progress: ~$sampling_lines/$num_samples sampling (~$pct%)"
        else
            pct=$((data_lines * 100 / num_warmup))
            echo "Progress: ~$data_lines/$num_warmup warmup (~$pct%)"
        fi
        echo ""

        # Current state (last 5 iterations)
        echo "Last 5 iterations:"
        tail -c 150000000 "$FILE" | grep -v "^#" | tail -5 | awk -F',' '{printf "  lp=%10.2f step=%.6f tree=%2d div=%d energy=%10.1f\n", $1, $3, $4, $6, $7}'
        echo ""

        # Divergence summary (recent only - full counts require full-file scan)
        echo "Divergence summary (recent):"
        divs_recent=$(tail -c 1500000000 "$FILE" | grep -v "^#" | tail -50 | cut -d',' -f6 | grep -c "^1" 2>/dev/null | tr -d '\n' || echo "0")
        echo "  Last 50 iterations: $divs_recent divergences"
        divs_recent_200=$(tail -c 5500000000 "$FILE" | grep -v "^#" | tail -200 | cut -d',' -f6 | grep -c "^1" 2>/dev/null | tr -d '\n' || echo "0")
        echo "  Last 200 iterations: $divs_recent_200 divergences"
        echo ""

        # Step size adaptation (warmup data is in the header region - read limited bytes)
        # For files with ~15MB lines, warmup of 300 iters = ~4.5GB, so we can only show
        # step sizes from the last N iterations, not the full warmup history
        echo "Recent step size history (last 200 iters, every 50):"
        tail -c 5500000000 "$FILE" | grep -v "^#" | tail -200 | awk -F',' '(NR%50)==0 {printf "  recent iter %3d: step=%.6f\n", NR, $3}' | head -10
        echo ""

        # Max treedepth (recent only)
        maxtree_recent=$(tail -c 5500000000 "$FILE" | grep -v "^#" | tail -200 | cut -d',' -f4 | grep -c "^10" 2>/dev/null | tr -d '\n' || echo "0")
        echo "Max treedepth hits (last 200): $maxtree_recent"
        ;;

    inits)
        echo "=== Initialization Comparison ==="
        echo ""

        # Key parameters to compare
        PARAMS="measure_sd tr_intercept_pop tr_sd_patient_intercept frac_intercept_pop init_intercept_pop"

        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                init_file=$(grep "^# init = " "$FILE" | head -1 | sed 's/.*= //' | tr -d ' ')
                echo "Chain $i init: $init_file"
                if [ -f "$init_file" ]; then
                    echo "  Key parameters:"
                    for param in $PARAMS; do
                        # Extract value using grep and sed (works for simple JSON)
                        val=$(grep -o "\"$param\"[[:space:]]*:[[:space:]]*[^,}]*" "$init_file" 2>/dev/null | head -1 | sed 's/.*://' | tr -d ' "')
                        if [ -n "$val" ]; then
                            printf "    %-30s = %s\n" "$param" "$val"
                        fi
                    done
                else
                    echo "  (init file not found or deleted)"
                fi
                echo ""
            fi
        done

        echo "Note: Init files in /tmp may be deleted after run completion."
        echo "To preserve inits, copy them before the run finishes."
        ;;

    cache)
        # Cache management - doesn't need RUN_ID
        subcommand="${CHAIN_NUM:-list}"
        case "$subcommand" in
            list)
                echo "=== Cached Diagnostics ==="
                echo "Cache directory: $CACHE_DIR"
                echo ""
                if [ -d "$CACHE_DIR" ]; then
                    for store_dir in "$CACHE_DIR"/*/; do
                        [ -d "$store_dir" ] || continue
                        store=$(basename "$store_dir")
                        echo "Store: $store"
                        for target_dir in "$store_dir"/*/; do
                            [ -d "$target_dir" ] || continue
                            target=$(basename "$target_dir")
                            for cache_file in "$target_dir"/*.cache; do
                                [ -f "$cache_file" ] || continue
                                run_id=$(basename "$cache_file" .cache)
                                # Get cache age
                                cache_mtime=$(stat -c %Y "$cache_file" 2>/dev/null || echo "0")
                                now=$(date +%s)
                                age_secs=$((now - cache_mtime))
                                if [ $age_secs -lt 60 ]; then
                                    age="${age_secs}s ago"
                                elif [ $age_secs -lt 3600 ]; then
                                    age="$((age_secs / 60))m ago"
                                else
                                    age="$((age_secs / 3600))h ago"
                                fi
                                echo "  $target / $run_id ($age)"
                            done
                        done
                    done
                else
                    echo "(no caches)"
                fi
                ;;
            clear)
                echo "Clearing all caches..."
                rm -rf "$CACHE_DIR"
                echo "Done."
                ;;
            *)
                echo "Unknown cache subcommand: $subcommand"
                echo "Usage: $0 cache [list|clear]"
                ;;
        esac
        exit 0
        ;;

    *)
        echo "Unknown command: $COMMAND"
        echo "Usage: $0 <output_dir> [command]"
        echo "Run '$0 --help' for available commands"
        exit 1
        ;;
esac
