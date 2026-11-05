# Genome-wide perturbation analysis

Snakemake workflow for genome-wide Perturb-seq screens (raw-count H5AD files).
The cells are split into chunks of perturbations, each chunk is run through
pertpy (nearest-control signatures) and Mixscape, and the cell labels are merged
afterwards. The merged data is then embedded and clustered once. Optionally,
a global DE compares all assigned cells against the controls.

## A note on the statistics

DE is done with a Wilcoxon rank-sum test on single cells. The reason is that
genome-wide screens usually have no biological replicates per perturbation: one
experiment, guides pooled, and often only a handful of cells per target. A
pseudobulk test (the better choice when you have replicates) needs several
independent samples per perturbation, and here we don't have them, so the cells
are the only units left to test on.

Be careful with this: cells from the same guide, batch or lane are not
independent, so the p-values are too optimistic. Use them to rank genes and
perturbations, and don't treat them as proof that a hit is real. If a dataset
does have replicates, use the optional pseudobulk tests instead.

Chunks are only there to split the compute. They are not replicates.

## Install

```bash
conda create -n xatlas --file workflow/envs/locks/linux-64.conda.txt
conda activate xatlas
python -m pip install --no-deps -r workflow/envs/locks/pip.txt
python -m unittest discover -s tests -p test_workflow_contracts.py -v
python tests/run_tiny_workflow.py --outdir examples/tiny
```

This is a Linux x86_64 environment (Snakemake 9.23.1, pertpy 1.1.1, Scanpy,
AnnData, Slurm plugin). I have not tested a fresh install or the GitHub Actions
yet. The last two commands run a small synthetic example through the whole
workflow.

## Run a screen

```bash
cp config/config.example.yaml config/my_screen.local.yaml
# edit dataset, paths, settings and results_dir
bash workflow/run_local.sh config/my_screen.local.yaml --cores 4 --dry-run
bash workflow/run_local.sh config/my_screen.local.yaml --cores 4
```

Each dataset needs `h5ad`, `pert_col` and `control_label`. `X` must be raw integer
counts, and cell and gene IDs must be unique. Perturbation labels come from
`obs[pert_col]` or from a `guide_calls_csv` (`cell_barcode,gene_target`, one row
per cell, so sort out multiple guides per cell beforehand). Write the dataset
version or accession into `data_id`.

Chunk composition, seeds and cell order change the PCA and the neighbours, so keep
them fixed if you want to reproduce a run. All chunks share one sampled set of
controls. A perturbation with more cells than the chunk size gets its own chunk.

On Slurm:

```bash
SLURM_ACCOUNT=your_account SLURM_PARTITION=your_partition \
  bash workflow/run_slurm.sh config/my_screen.local.yaml --dry-run
```

Keep the controller running in an allocation, and adjust memory and runtime with
Snakemake resource overrides.

## Outputs

```text
results_dir/
├── run_provenance.json     # settings, versions and file hashes of the run
└── DATASET/
    ├── inspection/         # checks of the input file
    ├── manifest.json, chunks/   # which perturbations are in which chunk, shared controls
    ├── chunk_runs/         # per-chunk labels and diagnostics
    ├── merged/             # selected cells
    ├── postprocess/        # pseudobulk profiles, UMAP/Leiden, cell-level DE
    └── global_de/          # optional DE of all assigned cells vs controls
```

Counts of the selected cells are summed per perturbation, normalised to 10,000
counts and log1p transformed. PCA spaces from different chunks cannot be
compared, so clustering is only done on the merged data in `postprocess/`.
The cluster table in `merged/` is empty if there is more than one chunk.

## Optional parts

Mixscale, R perturbation scoring, ORA, method comparison and the pseudobulk /
distribution tests (DESeq2, energy distance, Anderson-Darling) are optional and
switched on in the config. These are not validated on large atlas runs yet, and
the R environment is not locked. The full description of the pseudobulk and
distribution tests is in `docs/README_long_version_2026-11-05.md`.
