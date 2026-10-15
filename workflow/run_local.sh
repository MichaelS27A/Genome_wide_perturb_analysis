#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: bash workflow/run_local.sh CONFIG [Snakemake arguments]" >&2
  exit 2
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIGFILE="$(realpath "$1")"
shift
export PYTHONNOUSERSITE=1
exec snakemake --snakefile "$SCRIPT_DIR/Snakefile" \
  --directory "$REPO_ROOT" --configfile "$CONFIGFILE" \
  --cores 1 --rerun-incomplete "$@"
