#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

WORKDIR=$PWD \
SNAKEFILE=workflow/Snakefile \
CONFIGFILE=config/config.yaml \
SLURM_PARTITION=shortq \
SLURM_QOS=shortq \
SLURM_ACCOUNT=your_slurm_account \
TARGET_RULE=results/mixscape_pipeline/HCT116/chunk_summary.tsv \
JOBS=20 \
LOCAL_CORES=2 \
sbatch --partition=mediumq --qos=mediumq --account=your_slurm_account workflow/run_snakemake_master_longq_executor_slurm.sh
