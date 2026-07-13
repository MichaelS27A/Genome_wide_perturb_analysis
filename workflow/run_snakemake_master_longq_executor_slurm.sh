#!/usr/bin/env bash
#SBATCH --job-name=snakemake_master_slurmexec
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=lab_gsf
#SBATCH --output=snakemake_master_slurmexec_%j.out
#SBATCH --error=snakemake_master_slurmexec_%j.err

set -euo pipefail

# Master on current node (when run via sbatch) with child jobs submitted to SLURM.
# Behavior is aligned with workflow/run_interactive_snakemake.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
RAW_WORKDIR="${WORKDIR:-${SLURM_SUBMIT_DIR:-$REPO_ROOT}}"

# Accept submission from either repo root (.../code) or workflow dir (.../code/workflow).
if [[ -f "$RAW_WORKDIR/workflow/Snakefile" ]]; then
  WORKDIR="$RAW_WORKDIR"
elif [[ -f "$RAW_WORKDIR/Snakefile" && -d "$RAW_WORKDIR/rules" && -d "$RAW_WORKDIR/scripts" ]]; then
  WORKDIR="$(cd "$RAW_WORKDIR/.." && pwd)"
else
  WORKDIR="$REPO_ROOT"
fi
cd "$WORKDIR"

# ---- settings (aligned with interactive launcher) ----
SNAKEFILE="workflow/Snakefile"
CONFIGFILE="config/config.yaml"
TARGET_RULE="${TARGET_RULE:-all}"

# Keep this as the only easy override when submitting via sbatch:
#   sbatch --export=ALL,JOBS=10 workflow/run_snakemake_master_longq_executor_slurm.sh
JOBS="${JOBS:-50}"
LOCAL_CORES="${LOCAL_CORES:-1}"

# Child jobs default to shortq; master job queue is controlled by SBATCH lines above.
SLURM_PARTITION="${SLURM_PARTITION:-shortq}"
SLURM_QOS="${SLURM_QOS:-shortq}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lab_gsf}"

DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-8000}"
DEFAULT_RUNTIME_MIN="${DEFAULT_RUNTIME_MIN:-660}"
LATENCY_WAIT="${LATENCY_WAIT:-60}"
BUILD_CHUNK_MANIFESTS_MEM_MB="${BUILD_CHUNK_MANIFESTS_MEM_MB:-}"
# Default to mtime-only reruns for debug iterations; override when strict
# provenance-based reruns are desired.
RERUN_TRIGGERS="${RERUN_TRIGGERS:-mtime}"

CONDA_BASE="/nobackup/lab_gsf/mschoeber/miniconda3"
CONDA_ENV="snakemake"
# ------------------------------------------------------

if [[ ! -f "$SNAKEFILE" ]]; then
  echo "ERROR: missing $SNAKEFILE" >&2
  exit 1
fi
if [[ ! -f "$CONFIGFILE" ]]; then
  echo "ERROR: missing $CONFIGFILE" >&2
  exit 1
fi
if [[ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
  echo "ERROR: conda.sh not found at $CONDA_BASE/etc/profile.d/conda.sh" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

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
  --latency-wait "$LATENCY_WAIT"
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

echo "Running master; child jobs go to SLURM."
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"
