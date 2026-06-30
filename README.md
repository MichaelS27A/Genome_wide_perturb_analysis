# x_atlas perturbation workflow

this repository contains a snakemake workflow for large-scale perturb-seq analysis across multiple datasets.

goal

identify what genetic perturbations do at single-cell resolution, and prioritize perturbations that show reproducible effects across cell lines and cell types in public datasets.

implemented methods

- mixscape-based perturbation detection:
  - cells are grouped by perturbation and matched to control cells.
  - chunked execution is used to scale to large screens.
  - chunk-level outputs include perturbation effect summaries, effect vectors in pca space, and predicted perturbed cells.
- merged perturbation profiling:
  - chunk outputs are merged to dataset-level effect tables.
  - perturbations are clustered by effect-vector similarity.
- post-mixscape differential analysis:
  - selected perturbed cells are aggregated into perturbation-level pseudobulk profiles.
  - embedding and clustering are computed on perturbation-level profiles (pca, neighbors, umap, leiden).
  - differential expression is run per perturbation versus controls (wilcoxon), with ranked marker outputs.
- optional method hooks:
  - mixscale and perturbation score modules are included and can be enabled in the config.

datasets

- the default config is set up for multiple datasets and currently includes:
  - hct116 dual-guide perturb-seq data
  - hek293t dual-guide perturb-seq data (disabled by default in the template config)
- additional datasets can be added through `config/config.yaml` by defining dataset entries with input paths and perturbation/control labels.

repository layout

- `config/`: dataset registry and pipeline parameters.
- `workflow/`: snakemake rules, scripts, schemas, environments, and cluster launchers.
- `src/`: reusable library code.

quick start

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 4
```

notes

- large raw data and generated results are excluded from git.
- cluster launcher scripts under `workflow/` are provided for slurm-based runs and can be adapted for local policies.
