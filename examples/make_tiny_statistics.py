"""Small replicated count screen for exercising the optional statistical branches."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


def make_example(outdir, mode="biological_sample"):
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(19)
    genes = ["TARGET_A", "TARGET_B"] + [f"G{i}" for i in range(2, 80)]
    baseline = rng.uniform(8, 30, len(genes))
    rows, values = [], []
    for sample in range(4):
        sample_mean = baseline * rng.lognormal(0, 0.1, len(genes))
        for target, n_guides in (("Non-Targeting", 4), ("TARGET_A", 2), ("TARGET_B", 2)):
            for guide in range(n_guides):
                mean = sample_mean.copy()
                if target != "Non-Targeting":
                    mean[genes.index(target)] *= 0.1
                if target == "TARGET_A":
                    mean[2:10] *= 4
                for cell in range(20):
                    rows.append(dict(target=target, guide=f"{target}-{guide}",
                                     sample=f"donor-{sample}", batch=f"capture-{sample}"))
                    values.append(rng.poisson(mean))
    obs = pd.DataFrame(rows, index=[f"cell-{i}" for i in range(len(rows))])
    ad.AnnData(sparse.csr_matrix(np.asarray(values, dtype=np.int64)), obs=obs,
               var=pd.DataFrame(index=genes)).write_h5ad(outdir / "counts.h5ad")
    config = {
        "datasets": {"TINY_STATS": {
            "h5ad": str(outdir / "counts.h5ad"), "pert_col": "target",
            "control_label": "Non-Targeting", "data_id": "synthetic-replicated-statistics-v1",
            "statistics": {"guide_col": "guide", "sample_col": "sample", "batch_col": "batch",
                           "context_cols": [], "pseudobulk_mode": mode, "paired": True},
        }},
        "methods": {"mixscape": {"enabled": False}, "pseudobulk_de": {"enabled": True},
                    "energy_distance": {"enabled": True},
                    "anderson_darling_permutation": {"enabled": True}},
        "postprocess": {"enabled": False}, "audit": {"enabled": False},
        "results_dir": str(outdir / "results"),
        "conda_env": "envs/preprocessing.yaml",
        "statistics_conda_env": "envs/statistics.yaml",
        "pseudobulk_de": {
            "min_cells_per_pool": 5, "min_replicates": 3, "min_gene_total_count": 10,
            "min_gene_pools": 2, "fdr_alpha": 0.05, "ntc_null_groups": 1,
            "dispersion_fit": "mean",
            "random_seed": 0, "cells_per_read": 64, "max_count_nnz": 1000000,
            "max_model_elements": 1000000, "n_resamples": 3, "resample_fraction": 0.85,
            "guide_level": mode != "kolf_resample",
        },
        "distribution_tests": {
            "min_cells_per_group": 8, "max_control_cells": 40, "max_fit_cells": 64,
            "n_top_genes": 16, "pca_dims": 3, "n_permutations": 99,
            "random_seed": 0, "ntc_null_groups": 1, "normalize_target_sum": 10000,
            "max_dense_elements": 1000000, "max_distance_elements": 1000000,
        },
    }
    path = outdir / "config.yaml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    return path


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--mode", choices=["biological_sample", "guide_batch", "kolf_resample"],
                        default="biological_sample")
    args = parser.parse_args()
    print(make_example(args.outdir, args.mode))
