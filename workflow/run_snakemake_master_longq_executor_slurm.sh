#!/usr/bin/env bash
#SBATCH --job-name=snakemake_master_slurmexec
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --time=11:50:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --account=lab_gsf
#SBATCH --output=snakemake_master_slurmexec_%j.out
#SBATCH --error=snakemake_master_slurmexec_%j.err

set -euo pipefail

# Master on current node (when run via sbatch) with child jobs submitted to SLURM.
# Behavior is aligned with workflow/run_interactive_snakemake.sh.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${WORKDIR:-${SLURM_SUBMIT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}}"
cd "$WORKDIR"

# ---- settings (aligned with interactive launcher) ----
SNAKEFILE="${SNAKEFILE:-workflow/Snakefile}"
CONFIGFILE="${CONFIGFILE:-config/config.yaml}"
TARGET_RULE="${TARGET_RULE:-all}"

JOBS="${JOBS:-50}"
LOCAL_CORES="${LOCAL_CORES:-1}"

SLURM_PARTITION="${SLURM_PARTITION:-shortq}"
SLURM_QOS="${SLURM_QOS:-shortq}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-lab_gsf}"

DEFAULT_MEM_MB="${DEFAULT_MEM_MB:-8000}"
DEFAULT_RUNTIME_MIN="${DEFAULT_RUNTIME_MIN:-660}"
LATENCY_WAIT="${LATENCY_WAIT:-60}"

CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
CONDA_ENV="${CONDA_ENV:-snakemake}"
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
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

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
  --latency-wait "$LATENCY_WAIT"
  --rerun-incomplete
  --keep-going
)

if [[ -n "$SLURM_QOS" ]]; then
  CMD+=(--slurm-qos "$SLURM_QOS")
fi
if [[ "$TARGET_RULE" != "all" ]]; then
  CMD+=("$TARGET_RULE")
fi

echo "Running master; child jobs go to SLURM."
printf ' %q' "${CMD[@]}"
echo

"${CMD[@]}"
