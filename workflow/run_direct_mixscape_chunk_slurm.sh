#!/usr/bin/env bash
#SBATCH --job-name=mixscape_chunk_direct
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --time=11:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=180G
#SBATCH --account=lab_gsf
#SBATCH --output=mixscape_chunk_direct_%j.out
#SBATCH --error=mixscape_chunk_direct_%j.err

set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate "${PERTPY_CONDA_ENV:-pertpy_env}"
exec /nobackup/lab_gsf/mschoeber/data/X_Atlas/code/workflow/run_python_only_mixscape.sh chunk "${DATASET:-HCT116}" "${CHUNK:-0052}"
