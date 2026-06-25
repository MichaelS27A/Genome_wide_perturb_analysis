# X_Atlas Workflow Repo

Public workflow repository for perturbation-analysis pipelines in X_Atlas.

## Repository layout

- `config/`: workflow configuration and dataset registry.
- `src/`: reusable source modules (Python/R) that are not Snakemake rule scripts.
- `workflow/`: Snakemake entrypoint, rules, env specs, schemas, and executable workflow scripts.

## Quick start

1. Create/edit `config/config.yaml` and `config/samples.tsv` for your datasets.
2. Prepare environments referenced in `workflow/envs/`.
3. Run Snakemake from repo root:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 4
```

## Notes

- Large datasets and generated results are intentionally excluded from git.
- Cluster helper scripts under `workflow/` are templates and should be adapted for your scheduler/account.
