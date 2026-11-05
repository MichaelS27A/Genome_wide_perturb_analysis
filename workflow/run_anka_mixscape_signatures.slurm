#!/usr/bin/env bash
#SBATCH --job-name=anka_signatures
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --account=lab_gsf
#SBATCH --array=0-3
#SBATCH --output=/nobackup/lab_gsf/mschoeber/data/X_Atlas/code/logs/anka_signatures_%A_%a.out
#SBATCH --error=/nobackup/lab_gsf/mschoeber/data/X_Atlas/code/logs/anka_signatures_%A_%a.err
set -euo pipefail
workspace=/nobackup/lab_gsf/mschoeber/data
datasets=(K562_GWPS_RAW_FULL HEK293T HCT116 STEMCELL)
dataset="${datasets[${SLURM_ARRAY_TASK_ID:-0}]}"
python_bin="$workspace/X_Atlas/code/.snakemake/conda/658db681142ce08c6c4baccb86671aad_/bin/python"
export MPLCONFIGDIR=/tmp/anka_signature_matplotlib
export NUMBA_CACHE_DIR="/tmp/anka_numba_${SLURM_JOB_ID:-local}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMBA_NUM_THREADS=4
"$python_bin" -u "$workspace/X_Atlas/code/workflow/scripts/51_cluster_custom_mixscape_signatures.py" \
    --dataset "$dataset" --workspace "$workspace" \
    --outdir "$workspace/X_Atlas/results/anka_mixscape_signature_clustering_2026-10-19/$dataset"
