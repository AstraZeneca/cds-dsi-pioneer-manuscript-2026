#!/bin/bash
# Stan Chain Diagnostics Helper Script
# Usage: ./stan-diagnose-helper.sh <output_dir> [command]
#
# Commands:
#   status    - Quick health check of all chains (default)
#   stuck     - Find iterations with extreme lp__
#   warmup    - Show warmup progression
#   monitor   - Start background monitoring

set -e

OUTPUT_DIR="${1:-.}"
COMMAND="${2:-status}"

# Find the latest run
find_latest_run() {
    ls -t "$OUTPUT_DIR"/sf-ssm-log-space-*.csv 2>/dev/null | head -1 | grep -oP '\d{12}' | head -1
}

RUN_ID=$(find_latest_run)
if [ -z "$RUN_ID" ]; then
    echo "No Stan output files found in $OUTPUT_DIR"
    exit 1
fi

echo "Run ID: $RUN_ID"
echo "Directory: $OUTPUT_DIR"
echo "---"

# Stan CSV columns:
# 1=lp__, 2=accept_stat__, 3=stepsize__, 4=treedepth__, 5=n_leapfrog__, 6=divergent__, 7=energy__

case "$COMMAND" in
    status)
        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                lines=$(wc -l < "$FILE")
                lp=$(tail -1 "$FILE" | cut -d',' -f1)
                stepsize=$(tail -1 "$FILE" | cut -d',' -f3)
                treedepth=$(tail -1 "$FILE" | cut -d',' -f4)
                energy=$(tail -1 "$FILE" | cut -d',' -f7)
                # Cumulative divergences (column 6)
                divs=$(tail -n +49 "$FILE" | cut -d',' -f6 | grep -c "^1$" 2>/dev/null || echo "0")
                # Recent divergences (last 50 iterations)
                divs_recent=$(tail -50 "$FILE" | cut -d',' -f6 | grep -c "^1$" 2>/dev/null || echo "0")
                # Max treedepth hits (column 4, value = 10)
                maxtree=$(tail -n +49 "$FILE" | cut -d',' -f4 | grep -c "^10$" 2>/dev/null || echo "0")
                adapted=$(grep -c "Adaptation terminated" "$FILE" 2>/dev/null || echo "0")

                # Calculate progress (assuming 48 header + 500 warmup + 500 sampling = 1048)
                if [ "$adapted" -eq 1 ]; then
                    pct=$(( (lines - 548) * 100 / 500 ))
                    stage="sampling"
                else
                    pct=$(( (lines - 48) * 100 / 500 ))
                    stage="warmup"
                fi

                # Warn if problematic
                warn=""
                if [[ "$lp" == *"e+"* ]] && [[ "$lp" == "-"* ]]; then
                    warn=" [EXTREME LP!]"
                fi
                if [[ "$stepsize" == *"e-"* ]]; then
                    exp=$(echo "$stepsize" | grep -oP 'e-\d+' | grep -oP '\d+')
                    if [ "$exp" -gt 10 ]; then
                        warn="$warn [TINY STEPSIZE!]"
                    fi
                fi

                echo "Chain $i: $lines lines ($pct% $stage), lp=$lp, step=$stepsize, energy=$energy, divs=$divs(+$divs_recent), maxtree=$maxtree, adapted=$adapted$warn"
            else
                echo "Chain $i: FILE NOT FOUND"
            fi
        done
        ;;

    stuck)
        echo "Finding iterations with lp__ < -100000..."
        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                echo "=== Chain $i ==="
                awk -F',' 'NR>48 && $1 < -100000 {print NR-48": lp="$1" step="$3" tree="$4" div="$6" energy="$7}' "$FILE" | head -10
            fi
        done
        ;;

    warmup)
        echo "Warmup progression (every 50 iterations)..."
        for i in 1 2 3 4; do
            FILE=$(ls "$OUTPUT_DIR"/sf-ssm-log-space-${RUN_ID}-${i}-*.csv 2>/dev/null | head -1)
            if [ -f "$FILE" ]; then
                echo "=== Chain $i ==="
                awk -F',' 'NR>48 && (NR-48)%50==0 {printf "iter %d: lp=%.2f step=%.6f tree=%d div=%d energy=%.1f\n", NR-48, $1, $3, $4, $6, $7}' "$FILE" | head -15
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
                    lines=$(wc -l < "$FILE")
                    lp=$(tail -1 "$FILE" | cut -d',' -f1)
                    stepsize=$(tail -1 "$FILE" | cut -d',' -f3)
                    energy=$(tail -1 "$FILE" | cut -d',' -f7)
                    divs=$(tail -n +49 "$FILE" | cut -d',' -f6 | grep -c "^1$" 2>/dev/null || echo "0")
                    divs_recent=$(tail -50 "$FILE" | cut -d',' -f6 | grep -c "^1$" 2>/dev/null || echo "0")
                    maxtree=$(tail -n +49 "$FILE" | cut -d',' -f4 | grep -c "^10$" 2>/dev/null || echo "0")
                    adapted=$(grep -c "Adaptation terminated" "$FILE" 2>/dev/null || echo "0")
                    pct=$((($lines - 48) * 100 / 500))
                    echo "Chain $i: $lines/548 ($pct%) lp=$lp step=$stepsize energy=$energy divs=$divs(+$divs_recent) maxtree=$maxtree adapted=$adapted" >> "$LOG"
                fi
            done
            echo "Logged at $(date)"
            sleep 300
        done
        ;;

    *)
        echo "Unknown command: $COMMAND"
        echo "Usage: $0 <output_dir> [status|stuck|warmup|monitor]"
        exit 1
        ;;
esac
