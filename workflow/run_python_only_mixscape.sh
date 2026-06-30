#!/usr/bin/env bash
set -euo pipefail

# Direct Python-only Mixscape runner (no Snakemake).
#
# Modes:
# 1) Existing chunk:
#    workflow/run_python_only_mixscape.sh chunk HCT116 0052
#
# 2) One perturbation vs random controls:
#    workflow/run_python_only_mixscape.sh onepert HCT116 TP53
#
# Environment overrides:
#   CONTROL_N=10000        # controls in onepert mode
#   CHUNK_MAX_CONTROLS=10000   # cap controls in chunk mode; 0 disables capping
#   CHUNK_CONTROL_SEED=0       # RNG seed for chunk control downsampling
#   BATCH_SIZE=1024        # passed to 05_run_mixscape_chunk.py
#   PERT_COL=gene_target
#   CONTROL_LABEL=Non-Targeting
#   MIXSCAPE_CONDA_ENV=snakemake   # fallback env for conda run
#   PYTHON_BIN=/path/to/python     # force interpreter

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

MODE="${1:-chunk}"
DATASET="${2:-HCT116}"
ARG3="${3:-0052}"

PERT_COL="${PERT_COL:-gene_target}"
CONTROL_LABEL="${CONTROL_LABEL:-Non-Targeting}"
CONTROL_N="${CONTROL_N:-10000}"
CHUNK_MAX_CONTROLS="${CHUNK_MAX_CONTROLS:-10000}"
BATCH_SIZE="${BATCH_SIZE:-1024}"
SEED="${SEED:-0}"
CHUNK_CONTROL_SEED="${CHUNK_CONTROL_SEED:-$SEED}"
PCA_DIMS="${PCA_DIMS:-20}"
MIXSCAPE_CONDA_ENV="${MIXSCAPE_CONDA_ENV:-snakemake}"
PYTHON_BIN="${PYTHON_BIN:-python}"

PYTHON_CMD=()

check_python_modules() {
  local -a cmd=("$@")
  "${cmd[@]}" -c "import anndata, scanpy, pertpy" >/dev/null 2>&1
}

choose_python() {
  if check_python_modules "$PYTHON_BIN"; then
    PYTHON_CMD=("$PYTHON_BIN")
    return
  fi

  if command -v conda >/dev/null 2>&1; then
    local -a candidates=("$MIXSCAPE_CONDA_ENV" "preprocessing_ultra_minimal")
    local env
    for env in "${candidates[@]}"; do
      if conda run -n "$env" python -c "import anndata, scanpy, pertpy" >/dev/null 2>&1; then
        PYTHON_CMD=(conda run --no-capture-output -n "$env" python)
        echo "[env] using conda env: $env"
        return
      fi
    done
  fi

  echo "ERROR: could not find a Python with anndata/scanpy/pertpy." >&2
  echo "Try one of:" >&2
  echo "  MIXSCAPE_CONDA_ENV=snakemake bash workflow/run_python_only_mixscape.sh chunk HCT116 0052" >&2
  echo "  MIXSCAPE_CONDA_ENV=preprocessing_ultra_minimal bash workflow/run_python_only_mixscape.sh chunk HCT116 0052" >&2
  echo "  PYTHON_BIN=/path/to/python bash workflow/run_python_only_mixscape.sh chunk HCT116 0052" >&2
  exit 1
}

print_preflight() {
  echo "[preflight] python command:"
  printf ' %q' "${PYTHON_CMD[@]}"
  echo
  "${PYTHON_CMD[@]}" - <<'PY'
import sys
import anndata
import scanpy
import pertpy
print("[preflight] python:", sys.executable)
print("[preflight] anndata:", anndata.__version__)
print("[preflight] scanpy:", scanpy.__version__)
print("[preflight] pertpy:", pertpy.__version__)
PY
}

if [[ "$DATASET" != "HCT116" ]]; then
  echo "ERROR: only HCT116 defaults are wired in this helper right now." >&2
  echo "Set custom paths manually and run scripts/05_run_mixscape_chunk.py directly." >&2
  exit 1
fi

H5AD="../HCT116_filtered_dual_guide_cells.h5ad"
GUIDE_CSV="../HCT116_filtered_guide_calls_per_cell.csv.gz"
if [[ -d "../results/mixscape_pipeline" ]]; then
  RESULTS_ROOT="../results/mixscape_pipeline"
elif [[ -d "results/mixscape_pipeline" ]]; then
  RESULTS_ROOT="results/mixscape_pipeline"
else
  RESULTS_ROOT="../results/mixscape_pipeline"
fi
BASE_OUT="${RESULTS_ROOT}/${DATASET}/direct_python"

if [[ ! -f "$H5AD" ]]; then
  echo "ERROR: missing $H5AD" >&2
  exit 1
fi
if [[ "$MODE" == "onepert" && ! -f "$GUIDE_CSV" ]]; then
  echo "ERROR: missing $GUIDE_CSV" >&2
  exit 1
fi

choose_python
print_preflight

run_chunk() {
  local chunk="$1"
  local chunk_cells="${RESULTS_ROOT}/${DATASET}/chunks/chunk_${chunk}_cells.tsv.gz"
  local outdir="${BASE_OUT}/chunk_${chunk}"
  local chunk_cells_for_run="$chunk_cells"
  if [[ ! -f "$chunk_cells" ]]; then
    echo "ERROR: missing $chunk_cells" >&2
    exit 1
  fi
  mkdir -p "$outdir"

  # Cap controls to avoid scipy/CSR overflow/segfault paths in very large chunks.
  chunk_cells_for_run="$(
    "${PYTHON_CMD[@]}" - "$chunk_cells" "$outdir" "$chunk" "$CHUNK_MAX_CONTROLS" "$CHUNK_CONTROL_SEED" <<'PY'
import sys
from pathlib import Path
import pandas as pd

chunk_in = Path(sys.argv[1])
outdir = Path(sys.argv[2])
chunk_id = sys.argv[3]
max_controls = int(sys.argv[4])
seed = int(sys.argv[5])

if max_controls <= 0:
    print(chunk_in)
    raise SystemExit(0)

df = pd.read_csv(chunk_in, sep="\t", compression="infer")
if "is_control" not in df.columns:
    print(f"[chunk-cap] NOTE: no is_control column in {chunk_in}; using original chunk file.", file=sys.stderr)
    print(chunk_in)
    raise SystemExit(0)

is_ctrl = df["is_control"].astype(int) == 1
n_ctrl = int(is_ctrl.sum())
n_total = int(df.shape[0])
n_pert = n_total - n_ctrl
print(
    f"[chunk-cap] chunk={chunk_id} rows={n_total} controls={n_ctrl} perturbed={n_pert} max_controls={max_controls}",
    file=sys.stderr,
)

if n_ctrl <= max_controls:
    print(chunk_in)
    raise SystemExit(0)

outdir.mkdir(parents=True, exist_ok=True)
out_path = outdir / f"chunk_{chunk_id}_cells_ctrl{max_controls}.tsv.gz"
ctrl = df[is_ctrl].sample(n=max_controls, random_state=seed, replace=False)
pert = df[~is_ctrl]
out = pd.concat([ctrl, pert], axis=0, ignore_index=True)
out.to_csv(out_path, sep="\t", index=False, compression="gzip")
print(f"[chunk-cap] wrote {out_path} controls={len(ctrl)} perturbed={len(pert)}", file=sys.stderr)
print(out_path)
PY
  )"

  "${PYTHON_CMD[@]}" workflow/scripts/05_run_mixscape_chunk.py \
    --h5ad "$H5AD" \
    --chunk-cells "$chunk_cells_for_run" \
    --output-dir "$outdir" \
    --pert-col "$PERT_COL" \
    --control-label "$CONTROL_LABEL" \
    --pca-dims "$PCA_DIMS" \
    --chunk-id "$chunk" \
    --batch-size "$BATCH_SIZE"
}

run_one_pert() {
  local pert="$1"
  local safe_pert
  safe_pert="$(echo "$pert" | tr '/ ' '__')"
  local tmp_chunk="${BASE_OUT}/tmp_${safe_pert}_cells.tsv.gz"
  local outdir="${BASE_OUT}/onepert_${safe_pert}"
  mkdir -p "${BASE_OUT}" "$outdir"

  "${PYTHON_CMD[@]}" - <<PY
import pandas as pd, numpy as np
from pathlib import Path
rng = np.random.default_rng(${SEED})
csv_path = Path("${GUIDE_CSV}")
out_path = Path("${tmp_chunk}")
pert = "${pert}"
control_label = "${CONTROL_LABEL}"
control_n = int("${CONTROL_N}")

df = pd.read_csv(csv_path, sep="\\t", compression="infer", usecols=["cell_barcode","gene_target"], dtype="string")
df = df.dropna(subset=["cell_barcode","gene_target"]).copy()
df["cell_barcode"] = df["cell_barcode"].astype(str)
df["gene_target"] = df["gene_target"].astype(str)

ctrl = df[df["gene_target"] == control_label][["cell_barcode","gene_target"]].drop_duplicates("cell_barcode")
pert_df = df[df["gene_target"] == pert][["cell_barcode","gene_target"]].drop_duplicates("cell_barcode")

if pert_df.empty:
    raise SystemExit(f"No cells found for perturbation: {pert}")
if ctrl.empty:
    raise SystemExit(f"No control cells found for label: {control_label}")

if len(ctrl) > control_n:
    take = rng.choice(ctrl.index.to_numpy(), size=control_n, replace=False)
    ctrl = ctrl.loc[take]

ctrl = ctrl.assign(is_control=1)
pert_df = pert_df.assign(is_control=0)
out = pd.concat([ctrl, pert_df], axis=0, ignore_index=True)
out.to_csv(out_path, sep="\\t", index=False, compression="gzip")
print(f"Wrote {out_path} with n_control={len(ctrl)} n_pert={len(pert_df)}")
PY

  "${PYTHON_CMD[@]}" workflow/scripts/05_run_mixscape_chunk.py \
    --h5ad "$H5AD" \
    --chunk-cells "$tmp_chunk" \
    --output-dir "$outdir" \
    --pert-col "$PERT_COL" \
    --control-label "$CONTROL_LABEL" \
    --pca-dims "$PCA_DIMS" \
    --chunk-id "onepert_${safe_pert}" \
    --batch-size "$BATCH_SIZE"
}

case "$MODE" in
  check)
    echo "[ok] preflight passed"
    ;;
  chunk)
    run_chunk "$ARG3"
    ;;
  onepert)
    run_one_pert "$ARG3"
    ;;
  *)
    echo "ERROR: mode must be 'chunk' or 'onepert'" >&2
    exit 1
    ;;
esac
