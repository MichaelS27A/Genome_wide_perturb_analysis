#!/usr/bin/env bash
#SBATCH --job-name=snakemake_master_slurmexec
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --time=11:50:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=your_slurm_account
#SBATCH --output=snakemake_master_slurmexec_%j.out
#SBATCH --error=snakemake_master_slurmexec_%j.err

set -euo pipefail

# -----------------------------
# User-tunable settings
# -----------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer the sbatch submit directory; fall back to script directory.
WORKDIR="${WORKDIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"
SNAKEFILE="${SNAKEFILE:-workflow/Snakefile}"
TARGET_RULE="${TARGET_RULE:-all}"
JOBS="${JOBS:-400}"
LOCAL_CORES="${LOCAL_CORES:-4}"
LATENCY_WAIT="${LATENCY_WAIT:-120}"
RESTART_TIMES="${RESTART_TIMES:-1}"
SNAKEMAKE_CONDA_PREFIX="${SNAKEMAKE_CONDA_PREFIX:-$WORKDIR/.snakemake/conda}"
CONFIGFILE="${CONFIGFILE:-}"

# Cluster defaults for worker jobs submitted by Snakemake
SLURM_PARTITION="${SLURM_PARTITION:-shortq}"
SLURM_QOS="${SLURM_QOS:-$SLURM_PARTITION}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-your_slurm_account}"

# Conda setup
CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"
CONDA_ENV="snakemake"

LOG_DIR="${LOG_DIR:-$WORKDIR/snakemake_pipeline/.snakemake/slurm_logs}"
mkdir -p "$LOG_DIR"

echo "======================================================="
echo "Snakemake SLURM-executor master started on $(hostname) at $(date)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "WORKDIR=$WORKDIR"
echo "SNAKEFILE=$SNAKEFILE"
echo "TARGET_RULE=$TARGET_RULE"
echo "JOBS=$JOBS LOCAL_CORES=$LOCAL_CORES"
echo "SNAKEMAKE_CONDA_PREFIX=$SNAKEMAKE_CONDA_PREFIX"
if [[ -n "$CONFIGFILE" ]]; then
  echo "CONFIGFILE=$CONFIGFILE"
fi
echo "SLURM_PARTITION=$SLURM_PARTITION SLURM_QOS=${SLURM_QOS:-<none>}"
echo "======================================================="

cd "$WORKDIR"

# Keep threaded libs aligned with Slurm request.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

if [[ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]]; then
  source "$CONDA_BASE/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
else
  echo "ERROR: conda.sh not found at $CONDA_BASE/etc/profile.d/conda.sh"
  exit 1
fi

if ! command -v snakemake >/dev/null 2>&1; then
  echo "ERROR: snakemake not found in active environment ($CONDA_ENV)."
  exit 1
fi

echo "Snakemake version: $(snakemake --version)"

CONFIG_ARGS=()
if [[ -n "$CONFIGFILE" ]]; then
  CONFIG_ARGS+=(--configfile "$CONFIGFILE")
fi

echo "[step 0] Unlocking working directory"
snakemake -s "$SNAKEFILE" "${CONFIG_ARGS[@]}" --unlock || true

# Preferred cluster submit wrapper for classic --cluster mode.
CLUSTER_SUBMIT_SCRIPT="${CLUSTER_SUBMIT_SCRIPT:-$WORKDIR/workflow/cluster_submit.sh}"

if snakemake --help 2>/dev/null | grep -q -- '--executor'; then
  echo "Using Snakemake SLURM executor mode (--executor slurm)."
  SNAKEMAKE_ARGS=(
    -s "$SNAKEFILE"
    --use-conda
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX"
    --executor slurm
    --jobs "$JOBS"
    --cores 1
    --local-cores "$LOCAL_CORES"
    --default-resources
      "slurm_partition=$SLURM_PARTITION"
      "slurm_account=$SLURM_ACCOUNT"
      "cpus_per_task=1"
      "mem_mb=8000"
      "runtime=60"
    --latency-wait "$LATENCY_WAIT"
    --restart-times "$RESTART_TIMES"
    --rerun-triggers mtime
    --rerun-incomplete
    --keep-going
    "${CONFIG_ARGS[@]}"
  )
  if [[ -n "$SLURM_QOS" ]]; then
    SNAKEMAKE_ARGS+=(--slurm-qos "$SLURM_QOS")
  fi
  if [[ "$TARGET_RULE" != "all" ]]; then
    SNAKEMAKE_ARGS+=("$TARGET_RULE")
  fi
  snakemake "${SNAKEMAKE_ARGS[@]}"
else
  echo "Snakemake does not support --executor; falling back to --cluster sbatch mode."
  if [[ -x "$CLUSTER_SUBMIT_SCRIPT" ]]; then
    CLUSTER_CMD="bash $CLUSTER_SUBMIT_SCRIPT {resources.slurm_partition} {resources.mem_mb} {resources.runtime} {resources.gpu} {rule} {dependencies} {jobscript}"
    echo "Using cluster submit wrapper: $CLUSTER_SUBMIT_SCRIPT"
  else
    CLUSTER_CMD="sbatch --partition={resources.slurm_partition} --account=$SLURM_ACCOUNT --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --output=$LOG_DIR/%x_%j.out --error=$LOG_DIR/%x_%j.err"
    if [[ -n "$SLURM_QOS" ]]; then
      CLUSTER_CMD="sbatch --partition={resources.slurm_partition} --qos=$SLURM_QOS --account=$SLURM_ACCOUNT --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.runtime} --output=$LOG_DIR/%x_%j.out --error=$LOG_DIR/%x_%j.err"
    fi
    echo "cluster_submit.sh not found/executable; using inline sbatch command."
  fi
  SNAKEMAKE_ARGS=(
    -s "$SNAKEFILE"
    --use-conda
    --conda-prefix "$SNAKEMAKE_CONDA_PREFIX"
    --jobs "$JOBS"
    --cores 1
    --local-cores "$LOCAL_CORES"
    --cluster "$CLUSTER_CMD"
    --latency-wait "$LATENCY_WAIT"
    --restart-times "$RESTART_TIMES"
    --rerun-triggers mtime
    --rerun-incomplete
    --keep-going
    "${CONFIG_ARGS[@]}"
  )
  if [[ "$TARGET_RULE" != "all" ]]; then
    SNAKEMAKE_ARGS+=("$TARGET_RULE")
  fi
  snakemake "${SNAKEMAKE_ARGS[@]}"
fi

echo "======================================================="
echo "Snakemake SLURM-executor master finished at $(date)"
echo "======================================================="
