#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: SLURM_ACCOUNT=... SLURM_PARTITION=... bash workflow/run_slurm.sh CONFIG [Snakemake arguments]" >&2
  exit 2
fi
: "${SLURM_ACCOUNT:?Set SLURM_ACCOUNT to your allocation}"
: "${SLURM_PARTITION:?Set SLURM_PARTITION to your compute partition}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIGFILE="$(realpath "$1")"
shift
export PYTHONNOUSERSITE=1
exec snakemake --snakefile "$SCRIPT_DIR/Snakefile" --directory "$REPO_ROOT" \
  --configfile "$CONFIGFILE" --executor slurm --jobs "${JOBS:-20}" \
  --local-cores 1 --default-resources "slurm_account=$SLURM_ACCOUNT" \
  "slurm_partition=$SLURM_PARTITION" --latency-wait 60 --rerun-incomplete "$@"
