"""Shared config, defaults, and helper functions."""

from pathlib import Path

BASE_DIR = Path(workflow.basedir).resolve()
ROOT_DIR = BASE_DIR.parent

DEFAULT_CONFIG = {
    "datasets": {
        "HCT116": {
            "enabled": True,
            "h5ad": str(ROOT_DIR / "data" / "HCT116_filtered_dual_guide_cells.h5ad"),
            "guide_calls_csv": str(ROOT_DIR / "data" / "HCT116_filtered_guide_calls_per_cell.csv.gz"),
            "pert_col": "gene_target",
            "control_label": "Non-Targeting",
        },
        "HEK293T": {
            "enabled": True,
            "h5ad": str(ROOT_DIR / "data" / "HEK293T_filtered_dual_guide_cells.h5ad"),
            "guide_calls_csv": str(ROOT_DIR / "data" / "HEK293T_filtered_guide_calls_per_cell.csv.gz"),
            "pert_col": "gene_target",
            "control_label": "Non-Targeting",
        },
    },
    "chunking": {
        "perturbations_per_chunk": 128,
        "min_cells_per_perturbation": 30,
        "max_controls_per_chunk": 50000,
        "csv_read_chunk_size": 300000,
    },
    "mixscape": {
        "pca_dims": 20,
        "write_subset_h5ad": False,
        "use_hvg_for_pca": False,
        "normalize_target_sum": 10000,
        "logfc_threshold": 0.10,
        "pval_cutoff": 0.05,
    },
    "mixscale": {
        "use_hvg_for_pca": False,
        "normalize_target_sum": 10000,
        "logfc_threshold": 0.10,
        "pval_cutoff": 0.05,
        "min_de_genes": 5,
        "max_de_genes": 100,
        "batch_size": 0,
    },
    "clustering": {
        "n_clusters": 20,
    },
    "postprocess": {
        "min_selected_cells": 20,
        "max_control_cells_for_de": 50000,
        "n_top_de_genes": 200,
        "random_seed": 0,
    },
    "audit": {
        "enabled": True,
        "require_ready": True,
        "require_raw_source": True,
    },
    "methods": {
        "bootstrap_install": True,
        "mixscape": {"enabled": True},
        "mixscale": {
            "enabled": False,
            "min_cells_per_perturbation": 30,
            "max_perturbations": 0,
            "max_cells": 0,
            "random_seed": 0,
        },
        "ps": {
            "enabled": False,
            "min_cells_per_perturbation": 30,
            "max_perturbations": 0,
            "max_cells": 0,
            "random_seed": 0,
        },
        "comparison": {"enabled": False},
    },
    "results_dir": str(ROOT_DIR / "results" / "mixscape_pipeline"),
    "conda_env": "envs/preprocessing.yaml",
    "r_env": "envs/r_mix_methods.yaml",
    "rscript_bin": "Rscript",
}


def deep_merge(a, b):
    out = dict(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = v
    return out


CFG = deep_merge(DEFAULT_CONFIG, config)


def resolve_from_root(pathlike):
    p = Path(pathlike)
    if p.is_absolute():
        return p
    return (ROOT_DIR / p).resolve()


for _ds_name, _ds_cfg in CFG.get("datasets", {}).items():
    if _ds_cfg.get("h5ad"):
        _ds_cfg["h5ad"] = str(resolve_from_root(_ds_cfg["h5ad"]))
    if _ds_cfg.get("guide_calls_csv"):
        _ds_cfg["guide_calls_csv"] = str(resolve_from_root(_ds_cfg["guide_calls_csv"]))

if CFG.get("results_dir"):
    CFG["results_dir"] = str(resolve_from_root(CFG["results_dir"]))

RESULTS_DIR = Path(CFG["results_dir"])
DATASETS = [ds for ds, dcfg in CFG["datasets"].items() if dcfg.get("enabled", True)]

# Resolve env specs relative to workflow base dir so includes under rules/ do not
# accidentally shift them to rules/envs/...
if "conda_env" in CFG and CFG["conda_env"]:
    _conda_env = Path(CFG["conda_env"])
    if not _conda_env.is_absolute():
        CFG["conda_env"] = str((BASE_DIR / _conda_env).resolve())

if "r_env" in CFG and CFG["r_env"]:
    _r_env = Path(CFG["r_env"])
    if not _r_env.is_absolute():
        CFG["r_env"] = str((BASE_DIR / _r_env).resolve())

if not DATASETS:
    raise ValueError("No enabled datasets found in config.datasets")


def method_enabled(name):
    return bool(CFG.get("methods", {}).get(name, {}).get("enabled", False))


def audit_enabled():
    return bool(CFG.get("audit", {}).get("enabled", True))


def r_bootstrap_enabled():
    return bool(CFG.get("methods", {}).get("bootstrap_install", True))


def r_bootstrap_marker():
    return str(RESULTS_DIR / "_env" / "r_method_packages.ok")


MIXSCAPE_ENABLED = method_enabled("mixscape")
MIXSCALE_ENABLED = method_enabled("mixscale")
PS_ENABLED = method_enabled("ps")
COMPARISON_ENABLED = method_enabled("comparison")

if COMPARISON_ENABLED:
    missing = [m for m in ("mixscape", "mixscale", "ps") if not method_enabled(m)]
    if missing:
        raise ValueError(f"methods.comparison requires enabled methods: {missing}")


def all_targets():
    outs = []

    if audit_enabled():
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "inspection" / "ready.ok"), dataset=DATASETS))

    if PS_ENABLED:
        if r_bootstrap_enabled():
            outs.append(r_bootstrap_marker())

    if MIXSCAPE_ENABLED:
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "merge_meta.json"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_stats_merged.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_effects_pca_merged.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "perturbation_clusters.tsv"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "selected_perturbed_cells_merged.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "merged" / "mixscape_selection_summary.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "postprocess_meta.json"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_pseudobulk.h5ad"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_umap_leiden.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_differential_genes.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_long_table.tsv.gz"), dataset=DATASETS))
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "postprocess" / "perturbation_gene_long_table.tsv.gz"), dataset=DATASETS))

    if MIXSCALE_ENABLED:
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "mixscale" / "done.txt"), dataset=DATASETS))

    if PS_ENABLED:
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "ps" / "done.txt"), dataset=DATASETS))

    if COMPARISON_ENABLED:
        outs.extend(expand(str(RESULTS_DIR / "{dataset}" / "comparison" / "comparison_summary.tsv"), dataset=DATASETS))

    return outs


def maybe_audit_ready_input(wc):
    if audit_enabled():
        return str(RESULTS_DIR / wc.dataset / "inspection" / "ready.ok")
    return []


def get_chunk_ids(wildcards):
    ck = checkpoints.build_chunk_manifests.get(dataset=wildcards.dataset)
    chunk_dir = Path(ck.output.chunkdir)
    chunk_ids = sorted(glob_wildcards(str(chunk_dir / "chunk_{chunk}_cells.tsv.gz")).chunk)
    if not chunk_ids:
        raise ValueError(f"No chunk manifests found in {chunk_dir}")
    return chunk_ids


def chunk_cells_path(wildcards):
    ck = checkpoints.build_chunk_manifests.get(dataset=wildcards.dataset)
    return str(Path(ck.output.chunkdir) / f"chunk_{wildcards.chunk}_cells.tsv.gz")


def chunk_stats_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "perturbation_stats.tsv"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


def chunk_effects_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "perturbation_effects_pca.tsv"),
        dataset=wildcards.dataset,
        chunk=ids,
    )


def chunk_selected_inputs(wildcards):
    ids = get_chunk_ids(wildcards)
    return expand(
        str(RESULTS_DIR / "{dataset}" / "chunk_runs" / "chunk_{chunk}" / "selected_perturbed_cells.tsv.gz"),
        dataset=wildcards.dataset,
        chunk=ids,
    )
