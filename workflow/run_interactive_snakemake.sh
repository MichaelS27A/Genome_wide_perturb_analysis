#!/usr/bin/env bash
set -euo pipefail

# Minimal launcher:
# - Runs Snakemake master in this terminal
# - Submits child jobs to SLURM via --executor slurm

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

# ---- edit these if needed ----
SNAKEFILE="workflow/Snakefile"
CONFIGFILE="config/config.yaml"
TARGET_RULE="all"

JOBS=50
LOCAL_CORES=1

SLURM_PARTITION="shortq"
SLURM_QOS="shortq"
SLURM_ACCOUNT="lab_gsf"

DEFAULT_MEM_MB=8000
DEFAULT_RUNTIME_MIN=660
# ------------------------------

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
  --keep-going
)

if [[ -n "$SLURM_QOS" ]]; then
  CMD+=(--slurm-qos "$SLURM_QOS")
fi
if [[ "$TARGET_RULE" != "all" ]]; then
  CMD+=("$TARGET_RULE")
fi

echo "Running master locally; child jobs go to SLURM."
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"
