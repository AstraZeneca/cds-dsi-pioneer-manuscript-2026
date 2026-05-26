#!/bin/bash
# publication_targets.sh
#
# Wrapper for run_targets.sh with project set to "publication"
# Examples:
#   ./publication_targets.sh                 # run full pipeline
#   ./publication_targets.sh -c              # run with shortcut = TRUE
#   ./publication_targets.sh -m '"a","b"'    # run specific targets
#   ./publication_targets.sh -d -m '"a"'     # dry-run: print command only
#
# For full options, see run_targets.sh

# Get the directory containing this script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Call run_targets.sh with publication project name
# Pass all arguments through, adding "publication" at the end
exec "${SCRIPT_DIR}/run_targets.sh" "$@" publication
