"""Export raw-count pools for optional DESeq2; no cell-level expression selection."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).parent))
from workflow_support import axis_hash
from statistics_support import (
    aggregate, close_data, open_data, regular_membership, resample_membership,
    versions, write_json,
)

import anndata as ad
import numpy as np


def run(dataset, settings, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mode = dataset["statistics"]["pseudobulk_mode"]
    if mode == "biological_sample" and not dataset["statistics"].get("sample_col"):
        raise ValueError("Biological-sample mode requires an explicitly declared sample_col")
    if mode in {"guide_batch", "kolf_resample"} and not dataset["statistics"].get("batch_col"):
        raise ValueError("Guide/batch and KOLF modes require an explicitly declared batch_col")
    adata, handle, meta = open_data(dataset)
    try:
        if mode == "kolf_resample":
            membership, pools = resample_membership(
                meta, adata.n_obs, dataset["control_label"], settings["random_seed"],
                settings["n_resamples"], settings["resample_fraction"], settings["ntc_null_groups"],
            )
        else:
            membership, pools = regular_membership(meta, adata.n_obs)
        counts = aggregate(adata, handle, membership, settings["cells_per_read"],
                           settings["max_count_nnz"])
        pools.index = [f"P{i:08d}" for i in range(len(pools))]
        pools.index.name = "pool_id"
        pools["library_size"] = np.asarray(counts.sum(axis=1)).ravel()
        pools["eligible"] = (pools["n_cells"] >= settings["min_cells_per_pool"]) & (pools["library_size"] > 0)
        pools["exclusion_reason"] = np.where(pools["eligible"], "", "small_or_empty_pool")
        pools["mode"] = mode
        ad.AnnData(counts, obs=pools, var=adata.var.copy()).write_h5ad(outdir / "counts.h5ad")
        pools.to_csv(outdir / "samples.tsv.gz", sep="\t")
        pools[["n_cells", "library_size", "eligible", "exclusion_reason"]].to_csv(
            outdir / "pool_qc.tsv.gz", sep="\t")
        write_json(outdir / "aggregation_meta.json", {
            "dataset": dataset, "settings": settings, "mode": mode,
            "input_cells": adata.n_obs, "assigned_cells": len(meta),
            "n_pools": len(pools), "n_genes": adata.n_vars,
            "input_cell_axis_sha256": axis_hash(adata.obs_names),
            "gene_axis_sha256": axis_hash(adata.var_names),
            "membership_nnz": membership.nnz,
            "membership_sha256": axis_hash([
                axis_hash(membership.indptr), axis_hash(membership.indices), axis_hash(membership.data)]),
            "experimental_unit": {"biological_sample": "declared biological sample",
                                  "guide_batch": "guide/capture pool",
                                  "kolf_resample": "overlapping cell resample"}[mode],
            "package_versions": versions(),
        })
    finally:
        close_data(adata, handle)


if __name__ == "__main__":
    run(json.loads(snakemake.params.dataset), json.loads(snakemake.params.settings),
        snakemake.params.outdir)
