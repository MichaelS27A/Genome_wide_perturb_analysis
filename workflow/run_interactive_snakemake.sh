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
TARGET_RULE="${TARGET_RULE:-all}"

JOBS="${JOBS:-50}"
LOCAL_CORES="${LOCAL_CORES:-1}"

SLURM_PARTITION="${SLURM_PARTITION:-shortq}"
SLURM_QOS="${SLURM_QOS:-shortq}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lab_gsf}"

DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-8000}"
DEFAULT_RUNTIME_MIN="${DEFAULT_RUNTIME_MIN:-660}"
BUILD_CHUNK_MANIFESTS_MEM_MB="${BUILD_CHUNK_MANIFESTS_MEM_MB:-}"
# Default to mtime-only reruns for debug iterations; override when strict
# provenance-based reruns are desired.
RERUN_TRIGGERS="${RERUN_TRIGGERS:-mtime}"
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

# Avoid user-site Python packages (e.g., ~/.local snakemake) shadowing
# workflow/runtime dependencies inside rule jobs.
export PYTHONNOUSERSITE=1

echo "[step] unlock"
snakemake -s "$SNAKEFILE" --configfile "$CONFIGFILE" --unlock || true

RERUN_TRIGGERS_NORMALIZED="${RERUN_TRIGGERS//,/ }"
read -r -a RERUN_TRIGGER_ARGS <<< "$RERUN_TRIGGERS_NORMALIZED"
if [[ "${#RERUN_TRIGGER_ARGS[@]}" -eq 0 ]]; then
  RERUN_TRIGGER_ARGS=("mtime")
fi

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
  --rerun-triggers "${RERUN_TRIGGER_ARGS[@]}"
  --rerun-incomplete
  --keep-going
)

if [[ -n "$SLURM_QOS" ]]; then
  CMD+=(--slurm-qos "$SLURM_QOS")
fi
if [[ -n "$BUILD_CHUNK_MANIFESTS_MEM_MB" ]]; then
  CMD+=(--set-resources "build_chunk_manifests:mem_mb=${BUILD_CHUNK_MANIFESTS_MEM_MB}")
fi
if [[ "$TARGET_RULE" != "all" ]]; then
  CMD+=("$TARGET_RULE")
fi

echo "Running master locally; child jobs go to SLURM."
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"
