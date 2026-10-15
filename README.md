# Genome-wide perturbation analysis

Reusable Snakemake workflow for raw-count Perturb-seq H5AD datasets. It builds
seeded perturbation chunks, computes nearest-control expression signatures
with pertpy, classifies cells with Mixscape and merges selected-cell labels.
Expression profiles are embedded and clustered once after all chunks merge.
Optional global DE compares all eligible assigned cells against controls.

DE uses **cell-level Wilcoxon tests**. Aggregated expression profiles are for
visualization; they do not supply biological replicates for pseudobulk inference.

## Install and verify

The exact installed environment is exported for Linux x86_64. It includes
Snakemake 9.23.1, pertpy 1.1.1 pinned to an immutable Git commit, Scanpy,
AnnData and the Slurm executor plugin.

```bash
conda create -n xatlas --file workflow/envs/locks/linux-64.conda.txt
conda activate xatlas
python -m pip install --no-deps -r workflow/envs/locks/pip.txt
python -m unittest discover -s tests -p test_workflow_contracts.py -v
python tests/run_tiny_workflow.py --outdir examples/tiny
```

The Conda lock fixes artifact URLs and checksums; pip versions are fixed but
wheels are not hash-locked. The environment YAMLs offer version-pinned solves
for `--use-conda`, with less control over native-library builds. A fresh package
download/install was not tested in this sandbox.

The synthetic example runs the real two-chunk Mixscape, postprocessing, global
DE and merge workflow twice, checks knockdown direction, and compares numerical
outputs. Nine regression checks and both runs passed in the installed environment;
a standalone Mixscale chunk also passed. Fresh installation and GitHub Actions
have not been exercised here.

## Run a screen

```bash
cp config/config.example.yaml config/my_screen.local.yaml
# Edit the dataset, paths, settings and results_dir.
bash workflow/run_local.sh config/my_screen.local.yaml --cores 4 --dry-run
bash workflow/run_local.sh config/my_screen.local.yaml --cores 4
```

Use an activated environment. Add `--use-conda` to provision rule environments
instead. CLI configs are authoritative: unspecified default atlases are never
added. Relative input/result paths resolve from the repository root; environment
YAML paths resolve from `workflow/`. Launchers work from another directory when
given their path and the config's path.

Each dataset requires `h5ad`, `pert_col`, and `control_label`. H5AD `X` must
contain finite, nonnegative raw integer counts and unique cell/gene identifiers.
Labels can be in `obs[pert_col]` or an optional `guide_calls_csv` containing one
`cell_barcode,gene_target` row per cell. Resolve multiple guide assignments before
analysis. Exclusions and minimum-cell thresholds precede chunk creation.

Record an immutable input accession/version in `datasets.<name>.data_id`.
An optional `h5ad_sha256` is recorded but not verified automatically. Cell order,
chunk settings, seeds and environment matter: changing chunk composition changes
PCA bases and reference neighbors. A single sampled control set is shared across
chunks; chunks are compute partitions, not biological replicates. An individual
perturbation larger than the chunk cell budget occupies its own chunk without
dropping its cells.

For Slurm, set site settings explicitly and keep the controller running in a
suitable allocation:

```bash
SLURM_ACCOUNT=your_account SLURM_PARTITION=your_partition \
  bash workflow/run_slurm.sh config/my_screen.local.yaml --dry-run
```

Use Snakemake resource overrides to size memory/runtime for your screen. Older
site-specific launchers are retained. The recommended launchers preserve code,
parameter, input, environment and timestamp rerun checks, and never automatically
unlock another run.

## Outputs

```text
results_dir/
├── run_provenance.json            # resolved settings, source hashes, versions
└── DATASET/
    ├── inspection/                # structure report and configured-input checks
    ├── manifest.json, chunks/     # target partition and shared sampled controls
    ├── chunk_runs/chunk_*/        # labels, diagnostics and effective settings
    ├── merged/                    # selected cells and chunk diagnostics
    ├── postprocess/               # expression profiles, one UMAP/Leiden, cell DE
    └── global_de/                 # optional all-assigned-cell DE
```

The provenance record contains resolved settings, shared source hashes, input
paths/sizes/mtimes and package versions. Input sizes and mtimes are not content
checksums. Each classifier chunk records effective PCA/signature dimensions,
batch size and ordered cell/gene identifier hashes. Missing or corrupt chunk
outputs are errors; header-only tables can represent no eligible results.

Selected-cell raw counts are summed per perturbation, normalized to a fixed
10,000-count library size and log1p transformed. The merged profiles receive
one dataset-wide embedding. Fewer than three profiles or features yield missing
embedding values rather than artificial coordinates.

PCA bases fitted independently in chunks are incompatible. Merged effect tables
retain `chunk_id` and serve as diagnostics. Their cluster table is empty for
multiple chunks; use the postprocess embedding for dataset clustering. Earlier
results need an explicit rerun; this update does not validate historical results.

Mixscale, R perturbation scoring, ORA, method comparison and target-versus-random
clustering remain optional modules. Large atlas runs, optional branches and Slurm
need separate deployment validation. The R environment is not locked, and runtime
R package installation defaults to disabled.
