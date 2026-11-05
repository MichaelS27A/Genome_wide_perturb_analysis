# Genome-wide perturbation analysis

Reusable Snakemake workflow for raw-count Perturb-seq H5AD datasets. It builds
seeded perturbation chunks, computes nearest-control expression signatures
with pertpy, classifies cells with Mixscape and merges selected-cell labels.
Expression profiles are embedded and clustered once after all chunks merge.
Optional global DE compares all eligible assigned cells against controls.

The default DE analysis uses **cell-level Wilcoxon tests**. Aggregated expression
profiles are for visualization; they do not supply biological replicates for
the optional pseudobulk inference described below.

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

## Optional pseudobulk and distribution tests

The experimental statistics branches use raw counts and all assigned,
QC-passing cells independently of Mixscape. Enable `methods.pseudobulk_de`,
`methods.energy_distance`, or `methods.anderson_darling_permutation` explicitly.
The existing Wilcoxon branch remains available.

Declare experimental metadata for each enabled dataset. These grouping columns
come from H5AD `obs`; an optional assignment CSV overrides target labels and
selects cells but does not supply the sample or guide columns:

```yaml
statistics:
  guide_col: sgRNA
  sample_col: donor
  batch_col: capture_batch
  context_cols: [stimulation]
  pseudobulk_mode: biological_sample
  paired: true
```

`biological_sample` sums guides and captures within each declared biological
sample/target/context, then fits paired sample blocks or an explicitly unpaired
comparison. `guide_batch` retains guide/sample/capture pools and blocks matched
sample/capture units when paired. Capture labels and guides do not establish
biological replication. `kolf_resample` creates three overlapping 85% cell
resamples with matched same-stratum controls; these are cell resamples, and the
current reproduction mode tests pooled targets rather than individual guides.
Overlapping resamples can inflate significance; inspect their held-out control
diagnostics before interpreting DEG counts.

Guide agreement excludes the directly targeted gene. Target labels must match
expression gene IDs, or `statistics.gene_symbol_col` must name their mapping
column in `var`. Unmapped targets receive a missing correlation and a QC reason
instead of a potentially self-knockdown-driven agreement score.

DESeq2 uses explicit count filters, positive-count size factors, Wald contrasts,
Cook's filtering/refitting and independent filtering. BH adjustment is within
each complete gene/contrast family. `pseudobulk_de.dispersion_fit` is `parametric`
by default; select `mean` explicitly when appropriate. A failed parametric trend
does not automatically change the model. Insufficient matched samples are
recorded in QC; unexpected fitting failures stop the job. `ntc_null_groups`
requests distinct held-out pairs of NTC guides, excluded from their respective
control comparisons; setting it to zero explicitly disables calibration.

Distribution tests use a frozen, seeded fitting cohort per context/sample/batch
stratum, variance-selected genes, log1p library normalization and PCA. This is
a generic implementation, **not exact reproduction of atlas preprocessing**.
Energy distance uses the biased V-statistic estimator and label permutations;
BH adjustment spans perturbation/stratum tests within a context, with separate
target and NTC-null families. Anderson–Darling uses permutation calibration
across the same selected gene set, not the published extended `kSamples`
calibration. Its effect column describes a mean normalized-count ratio; the
AD statistic itself measures an unsigned distribution change.

`distribution_tests.n_permutations` defaults to 10,000. Monte Carlo p-values
cannot resolve below `1 / (n_permutations + 1)`; choose adequate resolution for
the tested family. Explicit count, dense-array and distance-matrix budgets stop
oversized jobs rather than silently changing features or dropping cells.

Outputs are in `pseudobulk_counts/`, `pseudobulk_de/` and `distribution_tests/`,
including full results, QC, effective settings, guide agreement where estimable,
and saved PCA bases. ORA/clustering integration and atlas-scale validation of
these new branches are pending.

Generate an executable small example:

```bash
python examples/make_tiny_statistics.py --outdir /tmp/tiny_statistics
bash workflow/run_local.sh /tmp/tiny_statistics/config.yaml --use-conda
```

The optional environment is `workflow/envs/statistics.yaml`. Its principal
package versions are pinned; a complete resolved environment lock and fresh
installation validation are pending. The underlying analysis can also be tested
directly inside that environment:

```bash
python -m unittest discover -s tests -p test_statistics.py -v
python tests/run_tiny_statistics.py --outdir /tmp/tiny_statistics_test --repeat
```

To extend the existing locked core environment, install the six optional
packages with `python -m pip install --no-deps -r workflow/envs/locks/statistics-pip.txt`
and run `python -m pip check`. CI is configured to exercise the statistics rules
in this extended environment as well as the existing core workflow.
