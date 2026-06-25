#!/usr/bin/env python3
"""Merge chunk-level Mixscape outputs into dataset-level summaries."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist


def load_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chunk-stats", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-effects", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-selected-cells", nargs="+", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--n-clusters", type=int, default=20)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    stats_df = pd.concat([load_table(p) for p in args.chunk_stats], axis=0, ignore_index=True)
    effects_df = pd.concat([load_table(p) for p in args.chunk_effects], axis=0, ignore_index=True)
    selected_df = pd.concat([load_table(p) for p in args.chunk_selected_cells], axis=0, ignore_index=True)

    stats_df = stats_df.drop_duplicates(subset=["perturbation"], keep="first").sort_values("perturbation")
    effects_df = effects_df.drop_duplicates(subset=["perturbation"], keep="first").sort_values("perturbation")

    stats_out = args.outdir / "perturbation_stats_merged.tsv.gz"
    effects_out = args.outdir / "perturbation_effects_pca_merged.tsv.gz"
    selected_out = args.outdir / "selected_perturbed_cells_merged.tsv.gz"
    selection_summary_out = args.outdir / "mixscape_selection_summary.tsv.gz"
    stats_df.to_csv(stats_out, sep="\t", index=False, compression="gzip")
    effects_df.to_csv(effects_out, sep="\t", index=False, compression="gzip")
    selected_df.to_csv(selected_out, sep="\t", index=False, compression="gzip")

    if selected_df.empty:
        sel_summary = pd.DataFrame(
            columns=["chunk_id", "perturbation", "n_selected_perturbed_cells", "mean_selected_cells_per_perturbation"]
        )
        avg_selected_cells = np.nan
    else:
        sel_counts = (
            selected_df.groupby(["chunk_id", "perturbation"], as_index=False)
            .agg(n_selected_perturbed_cells=("cell_barcode", "nunique"))
            .sort_values(["chunk_id", "perturbation"])
        )
        per_pert = (
            selected_df.groupby("perturbation", as_index=False)
            .agg(n_selected_perturbed_cells=("cell_barcode", "nunique"))
            .sort_values("perturbation")
        )
        avg_selected_cells = float(per_pert["n_selected_perturbed_cells"].mean()) if not per_pert.empty else np.nan
        sel_summary = sel_counts.merge(
            per_pert.rename(columns={"n_selected_perturbed_cells": "n_selected_total_for_perturbation"}),
            on="perturbation",
            how="left",
        )
        sel_summary["mean_selected_cells_per_perturbation"] = avg_selected_cells

    sel_summary.to_csv(selection_summary_out, sep="\t", index=False, compression="gzip")

    # Cluster perturbations by PCA effect vectors.
    pc_cols = [c for c in effects_df.columns if c.startswith("delta_pc")]
    cluster_out = args.outdir / "perturbation_clusters.tsv"

    if len(pc_cols) >= 2 and effects_df.shape[0] >= 3:
        mat = effects_df[pc_cols].to_numpy(dtype=float)
        d = pdist(mat, metric="euclidean")
        z = linkage(d, method="average")
        k = min(args.n_clusters, effects_df.shape[0] - 1)
        cl = fcluster(z, t=k, criterion="maxclust")
        cluster_df = pd.DataFrame(
            {
                "perturbation": effects_df["perturbation"].values,
                "cluster": cl,
                "effect_l2_vs_control": effects_df["effect_l2_vs_control"].values,
            }
        ).sort_values(["cluster", "perturbation"])
        cluster_df.to_csv(cluster_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["perturbation", "cluster", "effect_l2_vs_control"]).to_csv(
            cluster_out, sep="\t", index=False
        )

    meta = {
        "n_chunk_stats_files": len(args.chunk_stats),
        "n_chunk_effect_files": len(args.chunk_effects),
        "n_chunk_selected_cell_files": len(args.chunk_selected_cells),
        "n_perturbations_merged": int(effects_df.shape[0]),
        "mean_selected_cells_per_perturbation": avg_selected_cells,
        "outputs": {
            "stats": str(stats_out),
            "effects": str(effects_out),
            "clusters": str(cluster_out),
            "selected_perturbed_cells": str(selected_out),
            "selection_summary": str(selection_summary_out),
        },
    }
    (args.outdir / "merge_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"Wrote {args.outdir / 'merge_meta.json'}")


if __name__ == "__main__":
    main()
