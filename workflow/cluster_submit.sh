#!/usr/bin/env bash
# Cluster submission wrapper for Snakemake --cluster mode.
# Usage pattern from launcher:
#   bash cluster_submit.sh {resources.slurm_partition} {resources.mem_mb} {resources.runtime} {resources.gpu} {rule} {dependencies} {jobscript}

set -euo pipefail

PARTITION="${1:-tinyq}"
MEM_MB="${2:-4000}"
RUNTIME="${3:-60}"
GPU="${4:-0}"
RULE="${5:-rule}"
DEPENDENCIES="${6:-}"
JOBSCRIPT="${7:-}"

if [[ -z "$JOBSCRIPT" ]]; then
  echo "ERROR: Missing jobscript argument" >&2
  exit 1
fi

SLURM_QOS_DEFAULT="${SLURM_QOS_DEFAULT:-}"
SLURM_ACCOUNT_DEFAULT="${SLURM_ACCOUNT_DEFAULT:-your_slurm_account}"

LOGS_DIR="${SLURM_LOGS_DIR:-$(pwd)/snakemake_pipeline/.snakemake/slurm_logs}"
mkdir -p "$LOGS_DIR"

SBATCH_ARGS=(
  --partition="$PARTITION"
  --account="$SLURM_ACCOUNT_DEFAULT"
  --mem="${MEM_MB}M"
  --time="$RUNTIME"
  --cpus-per-task=1
  --job-name="smk.${RULE}"
  --output="${LOGS_DIR}/${RULE}_%j.out"
  --error="${LOGS_DIR}/${RULE}_%j.err"
)

if [[ -n "$SLURM_QOS_DEFAULT" ]]; then
  SBATCH_ARGS+=(--qos="$SLURM_QOS_DEFAULT")
fi

if [[ "$GPU" =~ ^[0-9]+$ ]] && [[ "$GPU" -gt 0 ]]; then
  SBATCH_ARGS+=(--gres="gpu:${GPU}")
fi

if [[ -n "$DEPENDENCIES" ]] && [[ "$DEPENDENCIES" != "None" ]]; then
  SBATCH_ARGS+=(--dependency="afterok:${DEPENDENCIES// /:}")
fi

sbatch "${SBATCH_ARGS[@]}" "$JOBSCRIPT" | awk '{print $4}'
