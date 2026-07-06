#!/usr/bin/env bash
set -euo pipefail

# Run exactly one mixscape chunk in this terminal, with child job submitted to SLURM.
#
# Usage:
#   workflow/run_one_chunk_interactive.sh [DATASET] [CHUNK]
# Example:
#   workflow/run_one_chunk_interactive.sh HCT116 0052

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

DATASET="${1:-HCT116}"
CHUNK="${2:-0052}"
TARGET="results/mixscape_pipeline/${DATASET}/chunk_runs/chunk_${CHUNK}/done.txt"

SNAKEFILE="workflow/Snakefile"
CONFIGFILE="config/config.yaml"

LOCAL_CORES="${LOCAL_CORES:-1}"
JOBS="${JOBS:-1}"

SLURM_PARTITION="${SLURM_PARTITION:-shortq}"
SLURM_QOS="${SLURM_QOS:-shortq}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lab_gsf}"

DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-8000}"
DEFAULT_RUNTIME_MIN="${DEFAULT_RUNTIME_MIN:-660}"

if [[ ! -f "$SNAKEFILE" ]]; then
  echo "ERROR: missing $SNAKEFILE" >&2
  exit 1
fi
if [[ ! -f "$CONFIGFILE" ]]; then
  echo "ERROR: missing $CONFIGFILE" >&2
  exit 1
fi
if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake not found on PATH" >&2
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "ERROR: sbatch not found on PATH" >&2
  exit 1
fi

# Avoid user-site Python packages (e.g., ~/.local snakemake) shadowing
# workflow/runtime dependencies inside rule jobs.
export PYTHONNOUSERSITE=1

echo "[step] unlock"
snakemake -s "$SNAKEFILE" --configfile "$CONFIGFILE" --unlock || true

CMD=(
  snakemake
  -s "$SNAKEFILE"
  --configfile "$CONFIGFILE"
  --use-conda
  --executor slurm
  --jobs "$JOBS"
  --local-cores "$LOCAL_CORES"
  --default-resources
  "slurm_partition=$SLURM_PARTITION"
  "slurm_account=$SLURM_ACCOUNT"
  "mem_mb=$DEFAULT_MEM_MB"
  "runtime=$DEFAULT_RUNTIME_MIN"
  "cpus_per_task=1"
  --latency-wait 60
  --rerun-incomplete
  --forcerun build_chunk_manifests
  --keep-going
  "$TARGET"
)

if [[ -n "$SLURM_QOS" ]]; then
  CMD+=(--slurm-qos "$SLURM_QOS")
fi

echo "Running one chunk target:"
echo "  dataset=$DATASET chunk=$CHUNK"
echo "  target=$TARGET"
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"
