#!/usr/bin/env python3
"""Post-Mixscape perturbation analysis: pseudobulk, UMAP/Leiden, and DE."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _to_dense_row(x):
    if sparse.issparse(x):
        return np.asarray(x).ravel()
    return np.asarray(x).ravel()


def extract_rank_genes_groups(adata: sc.AnnData, groups: list[str], n_top: int) -> pd.DataFrame:
    rg = adata.uns["rank_genes_groups"]
    rows = []
    for g in groups:
        names = rg["names"][g][:n_top]
        scores = rg["scores"][g][:n_top]
        logfc = rg.get("logfoldchanges", None)
        pvals = rg.get("pvals", None)
        pvals_adj = rg.get("pvals_adj", None)
        pts = rg.get("pts", None)
        pts_rest = rg.get("pts_rest", None)

        for i, gene in enumerate(names, start=1):
            rec = {
                "perturbation": g,
                "rank": i,
                "gene": str(gene),
                "score": float(scores[i - 1]) if np.isfinite(scores[i - 1]) else np.nan,
                "logfoldchange": np.nan,
                "pval": np.nan,
                "pval_adj": np.nan,
                "pct_nz_group": np.nan,
                "pct_nz_reference": np.nan,
            }
            if logfc is not None:
                v = logfc[g][i - 1]
                rec["logfoldchange"] = float(v) if np.isfinite(v) else np.nan
            if pvals is not None:
                v = pvals[g][i - 1]
                rec["pval"] = float(v) if np.isfinite(v) else np.nan
            if pvals_adj is not None:
                v = pvals_adj[g][i - 1]
                rec["pval_adj"] = float(v) if np.isfinite(v) else np.nan
            if pts is not None and g in pts.columns and str(gene) in pts.index:
                rec["pct_nz_group"] = float(pts.loc[str(gene), g])
            if pts_rest is not None and g in pts_rest.columns and str(gene) in pts_rest.index:
                rec["pct_nz_reference"] = float(pts_rest.loc[str(gene), g])
            rows.append(rec)

    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--selected-cells", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--min-selected-cells", type=int, default=20)
    ap.add_argument("--max-control-cells", type=int, default=50000)
    ap.add_argument("--n-top-de-genes", type=int, default=200)
    ap.add_argument("--random-seed", type=int, default=0)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.random_seed)

    selected = pd.read_csv(args.selected_cells, sep="\t", compression="infer", dtype="string")
    if selected.empty:
        raise RuntimeError("selected-cells table is empty")

    selected = selected.dropna(subset=["cell_barcode", "perturbation"])
    selected["cell_barcode"] = selected["cell_barcode"].astype(str)
    selected["perturbation"] = selected["perturbation"].astype(str)

    pert_counts = (
        selected.groupby("perturbation", as_index=False)
        .agg(n_selected_cells=("cell_barcode", "nunique"))
        .sort_values("n_selected_cells", ascending=False)
    )
    keep_perts = pert_counts.loc[pert_counts["n_selected_cells"] >= args.min_selected_cells, "perturbation"]
    selected = selected[selected["perturbation"].isin(set(keep_perts))].copy()

    if selected.empty:
        raise RuntimeError(
            f"No perturbations with >= {args.min_selected_cells} selected cells after Mixscape filtering"
        )

    selected_map = selected.drop_duplicates(subset=["cell_barcode"]).set_index("cell_barcode")["perturbation"]

    adata_backed = sc.read_h5ad(args.h5ad, backed="r")
    obs = adata_backed.obs[[args.pert_col]].copy()
    obs.index = obs.index.astype(str)
    obs[args.pert_col] = obs[args.pert_col].astype(str)

    control_barcodes = obs.index[obs[args.pert_col] == args.control_label].to_numpy(dtype=str)
    if control_barcodes.size == 0:
        adata_backed.file.close()
        raise RuntimeError(f"No control cells found for label '{args.control_label}'")

    if args.max_control_cells > 0 and control_barcodes.size > args.max_control_cells:
        control_barcodes = rng.choice(control_barcodes, size=args.max_control_cells, replace=False)

    keep_barcodes = pd.Index(selected_map.index.tolist() + control_barcodes.tolist()).unique()
    obs_names = pd.Index(adata_backed.obs_names.astype(str))
    keep_mask = obs_names.isin(keep_barcodes)

    if int(keep_mask.sum()) == 0:
        adata_backed.file.close()
        raise RuntimeError("No overlap between selected/control barcodes and adata.obs_names")

    adata = adata_backed[keep_mask].to_memory()
    adata_backed.file.close()

    adata.obs_names = adata.obs_names.astype(str)
    adata.obs = adata.obs.copy()

    adata.obs["selected_perturbation"] = adata.obs_names.map(selected_map).astype("string")
    adata.obs["de_group"] = np.where(
        adata.obs["selected_perturbation"].notna(),
        adata.obs["selected_perturbation"].astype(str),
        args.control_label,
    )

    # Build pseudobulk matrix from selected perturbed cells only.
    adata_sel = adata[adata.obs["selected_perturbation"].notna()].copy()
    X = adata_sel.layers["counts"] if "counts" in adata_sel.layers else adata_sel.X

    perts = sorted(adata_sel.obs["selected_perturbation"].astype(str).unique().tolist())
    pidx = {p: i for i, p in enumerate(perts)}

    pb_rows = []
    pb_meta = []
    for p in perts:
        mask = adata_sel.obs["selected_perturbation"].astype(str).values == p
        if int(mask.sum()) == 0:
            continue
        row = X[mask].sum(axis=0)
        pb_rows.append(_to_dense_row(row))
        pb_meta.append({"perturbation": p, "n_selected_cells": int(mask.sum())})

    pb_X = np.vstack(pb_rows)
    pb_obs = pd.DataFrame(pb_meta).set_index("perturbation")
    pb = sc.AnnData(X=pb_X, obs=pb_obs, var=adata_sel.var.copy())

    sc.pp.normalize_total(pb)
    sc.pp.log1p(pb)

    if pb.n_obs >= 3:
        sc.pp.highly_variable_genes(pb, n_top_genes=min(2000, pb.n_vars), subset=False)
        sc.pp.pca(pb)
        sc.pp.neighbors(pb, n_neighbors=min(15, pb.n_obs - 1))
        sc.tl.umap(pb)
        sc.tl.leiden(pb, key_added="leiden", flavor="igraph", n_iterations=2)
    else:
        pb.obsm["X_umap"] = np.zeros((pb.n_obs, 2), dtype=float)
        pb.obs["leiden"] = "0"

    umap_df = pb.obs.copy().reset_index().rename(columns={"index": "perturbation"})
    umap_df["umap_1"] = pb.obsm["X_umap"][:, 0]
    umap_df["umap_2"] = pb.obsm["X_umap"][:, 1]

    # Differential expression for each perturbation vs controls.
    adata_de = adata.copy()
    sc.pp.normalize_total(adata_de)
    sc.pp.log1p(adata_de)

    de_counts = adata_de.obs["de_group"].value_counts()
    valid_groups = [
        g for g in de_counts.index.tolist() if g != args.control_label and int(de_counts[g]) >= args.min_selected_cells
    ]

    if valid_groups:
        sc.tl.rank_genes_groups(
            adata_de,
            groupby="de_group",
            groups=valid_groups,
            reference=args.control_label,
            method="wilcoxon",
            pts=True,
        )
        de_df = extract_rank_genes_groups(adata_de, valid_groups, args.n_top_de_genes)
    else:
        de_df = pd.DataFrame(
            columns=[
                "perturbation",
                "rank",
                "gene",
                "score",
                "logfoldchange",
                "pval",
                "pval_adj",
                "pct_nz_group",
                "pct_nz_reference",
            ]
        )

    pert_summary = (
        selected.groupby(["chunk_id", "perturbation"], as_index=False)
        .agg(n_selected_cells=("cell_barcode", "nunique"))
        .sort_values(["chunk_id", "perturbation"])
    )

    long_table = umap_df.merge(
        pert_summary.groupby("perturbation", as_index=False)["n_selected_cells"].sum(),
        on="perturbation",
        how="left",
        suffixes=("", "_from_chunks"),
    )

    long_table.to_csv(args.outdir / "perturbation_long_table.tsv.gz", sep="\t", index=False, compression="gzip")
    gene_long = de_df.merge(
        umap_df[["perturbation", "leiden", "umap_1", "umap_2", "n_selected_cells"]],
        on="perturbation",
        how="left",
    )
    gene_long.to_csv(
        args.outdir / "perturbation_gene_long_table.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pert_summary.to_csv(args.outdir / "perturbation_chunk_selection_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    umap_df.to_csv(args.outdir / "perturbation_umap_leiden.tsv.gz", sep="\t", index=False, compression="gzip")
    de_df.to_csv(args.outdir / "perturbation_differential_genes.tsv.gz", sep="\t", index=False, compression="gzip")
    pb.write_h5ad(args.outdir / "perturbation_pseudobulk.h5ad", compression="gzip")

    meta = {
        "h5ad": str(args.h5ad),
        "selected_cells": str(args.selected_cells),
        "n_selected_cells_total": int(selected["cell_barcode"].nunique()),
        "n_perturbations_selected": int(selected["perturbation"].nunique()),
        "n_perturbations_pseudobulk": int(pb.n_obs),
        "n_controls_used_for_de": int((adata.obs["de_group"] == args.control_label).sum()),
        "n_de_rows": int(de_df.shape[0]),
        "min_selected_cells": int(args.min_selected_cells),
        "n_top_de_genes": int(args.n_top_de_genes),
    }
    (args.outdir / "postprocess_meta.json").write_text(json.dumps(meta, indent=2))

    print(f"Wrote {args.outdir / 'postprocess_meta.json'}")


if __name__ == "__main__":
    main()
