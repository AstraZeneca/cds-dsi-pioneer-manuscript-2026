#!/usr/bin/env bash
# Domino App entry point — serves the pre-rendered Pioneer Quarto website.
#
# The _site directory is written to /mnt/artifacts by pioneer_render.sh
# and served here without any build step at startup.

set -euo pipefail

site_dir="${DOMINO_ARTIFACTS_DIR}/${DOMINO_STARTING_USERNAME}/pioneer/website/_site"

if [ ! -d "${site_dir}" ]; then
    echo "Error: site directory not found at ${site_dir}" >&2
    echo "Run pioneer_render.sh first to build the site." >&2
    exit 1
fi

echo "Serving ${site_dir} on port 8888"
npx serve "${site_dir}" --listen 8888 --no-clipboard
