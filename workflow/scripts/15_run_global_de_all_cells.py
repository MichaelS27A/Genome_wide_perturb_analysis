#!/usr/bin/env python3
"""Run global DE per perturbation vs control using all assigned cells (optionally capped)."""

import argparse
import json
import os
from pathlib import Path

# Apply conservative defaults before importing numerical libraries.
for _env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_env, "1")

import numpy as np
import pandas as pd
import scanpy as sc


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


def sample_group_barcodes(
    obs: pd.DataFrame,
    pert_col: str,
    control_label: str,
    valid_groups: list[str],
    max_control_cells: int,
    max_cells_per_pert: int,
    seed: int,
) -> tuple[pd.Index, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    selected: list[str] = []

    ctrl = obs.index[obs[pert_col] == control_label].to_numpy(dtype=str)
    if max_control_cells > 0 and ctrl.size > max_control_cells:
        ctrl = rng.choice(ctrl, size=max_control_cells, replace=False)
    selected.extend(ctrl.tolist())

    for p in valid_groups:
        idx = obs.index[obs[pert_col] == p].to_numpy(dtype=str)
        if max_cells_per_pert > 0 and idx.size > max_cells_per_pert:
            idx = rng.choice(idx, size=max_cells_per_pert, replace=False)
        selected.extend(idx.tolist())

    keep = pd.Index(selected).unique()
    counts = (
        obs.loc[obs.index.intersection(keep), [pert_col]]
        .groupby(pert_col, as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
        .sort_values("n_cells", ascending=False)
        .reset_index(drop=True)
    )
    return keep, counts


def run_analysis(args: argparse.Namespace) -> None:
    args.h5ad = Path(args.h5ad)
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata_backed = sc.read_h5ad(args.h5ad, backed="r")
    obs = adata_backed.obs[[args.pert_col]].copy()
    obs.index = obs.index.astype(str)
    obs[args.pert_col] = obs[args.pert_col].astype(str)

    if args.control_label not in set(obs[args.pert_col].tolist()):
        adata_backed.file.close()
        raise RuntimeError(f"control label '{args.control_label}' not found in '{args.pert_col}'")

    counts_raw = (
        obs.groupby(args.pert_col, as_index=False)
        .size()
        .rename(columns={"size": "n_cells_raw"})
        .sort_values("n_cells_raw", ascending=False)
        .reset_index(drop=True)
    )
    valid_groups = [
        p
        for p, n in zip(counts_raw[args.pert_col].tolist(), counts_raw["n_cells_raw"].tolist())
        if p != args.control_label and int(n) >= int(args.min_cells_per_perturbation)
    ]
    if not valid_groups:
        adata_backed.file.close()
        empty = pd.DataFrame(
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
        empty.to_csv(args.outdir / "perturbation_differential_genes.tsv.gz", sep="\t", index=False, compression="gzip")
        counts_raw.to_csv(args.outdir / "group_cell_counts.tsv.gz", sep="\t", index=False, compression="gzip")
        meta = {
            "h5ad": str(args.h5ad),
            "pert_col": args.pert_col,
            "control_label": args.control_label,
            "n_cells_total": int(obs.shape[0]),
            "n_perturbations_valid": 0,
            "note": "No perturbations passed min_cells_per_perturbation threshold.",
        }
        (args.outdir / "de_meta.json").write_text(json.dumps(meta, indent=2))
        (args.outdir / "done.txt").write_text("ok\n")
        return

    keep_barcodes, counts_used = sample_group_barcodes(
        obs=obs,
        pert_col=args.pert_col,
        control_label=args.control_label,
        valid_groups=valid_groups,
        max_control_cells=int(args.max_control_cells),
        max_cells_per_pert=int(args.max_cells_per_perturbation),
        seed=int(args.random_seed),
    )

    obs_names = pd.Index(adata_backed.obs_names.astype(str))
    keep_mask = obs_names.isin(keep_barcodes)
    if int(keep_mask.sum()) == 0:
        adata_backed.file.close()
        raise RuntimeError("No overlap between selected barcodes and adata.obs_names")

    adata = adata_backed[keep_mask].to_memory()
    adata_backed.file.close()
    adata.obs_names = adata.obs_names.astype(str)
    adata.obs = adata.obs.copy()
    adata.obs[args.pert_col] = adata.obs[args.pert_col].astype(str)

    # Keep only groups present after optional capping/subsampling.
    counts_after = adata.obs[args.pert_col].value_counts()
    valid_after = [
        p
        for p in valid_groups
        if p in counts_after.index and int(counts_after[p]) >= int(args.min_cells_per_perturbation)
    ]
    if not valid_after:
        empty = pd.DataFrame(
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
        empty.to_csv(args.outdir / "perturbation_differential_genes.tsv.gz", sep="\t", index=False, compression="gzip")
        counts_used.to_csv(args.outdir / "group_cell_counts.tsv.gz", sep="\t", index=False, compression="gzip")
        meta = {
            "h5ad": str(args.h5ad),
            "pert_col": args.pert_col,
            "control_label": args.control_label,
            "n_cells_total": int(adata.n_obs),
            "n_perturbations_valid": 0,
            "note": "No perturbations remained above threshold after optional capping.",
        }
        (args.outdir / "de_meta.json").write_text(json.dumps(meta, indent=2))
        (args.outdir / "done.txt").write_text("ok\n")
        return

    if float(args.normalize_target_sum) > 0:
        sc.pp.normalize_total(adata, target_sum=float(args.normalize_target_sum))
    if bool(args.log1p):
        sc.pp.log1p(adata)

    n_top = int(args.n_top_de_genes)
    if n_top <= 0:
        n_top = int(adata.n_vars)
    n_top = min(n_top, int(adata.n_vars))

    sc.tl.rank_genes_groups(
        adata,
        groupby=args.pert_col,
        groups=valid_after,
        reference=args.control_label,
        method=str(args.method),
        n_genes=n_top,
        pts=True,
    )
    de_df = extract_rank_genes_groups(adata, groups=valid_after, n_top=n_top)
    de_df.to_csv(args.outdir / "perturbation_differential_genes.tsv.gz", sep="\t", index=False, compression="gzip")
    counts_used.to_csv(args.outdir / "group_cell_counts.tsv.gz", sep="\t", index=False, compression="gzip")

    meta = {
        "h5ad": str(args.h5ad),
        "pert_col": args.pert_col,
        "control_label": args.control_label,
        "method": args.method,
        "normalize_target_sum": float(args.normalize_target_sum),
        "log1p": bool(args.log1p),
        "n_top_de_genes": int(n_top),
        "min_cells_per_perturbation": int(args.min_cells_per_perturbation),
        "max_control_cells": int(args.max_control_cells),
        "max_cells_per_perturbation": int(args.max_cells_per_perturbation),
        "random_seed": int(args.random_seed),
        "n_cells_used": int(adata.n_obs),
        "n_genes_used": int(adata.n_vars),
        "n_perturbations_valid": int(len(valid_after)),
        "n_de_rows": int(de_df.shape[0]),
    }
    (args.outdir / "de_meta.json").write_text(json.dumps(meta, indent=2))
    (args.outdir / "done.txt").write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--method", type=str, default="wilcoxon")
    ap.add_argument("--normalize-target-sum", type=float, default=10000)
    ap.add_argument("--log1p", action="store_true")
    ap.add_argument("--n-top-de-genes", type=int, default=2000)
    ap.add_argument("--min-cells-per-perturbation", type=int, default=30)
    ap.add_argument("--max-control-cells", type=int, default=0)
    ap.add_argument("--max-cells-per-perturbation", type=int, default=0)
    ap.add_argument("--random-seed", type=int, default=0)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        h5ad=Path(str(snk.input.h5ad)),
        outdir=Path(str(snk.params.outdir)),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        method=str(snk.params.method),
        normalize_target_sum=float(snk.params.normalize_target_sum),
        log1p=bool(snk.params.log1p),
        n_top_de_genes=int(snk.params.n_top_de_genes),
        min_cells_per_perturbation=int(snk.params.min_cells),
        max_control_cells=int(snk.params.max_control_cells),
        max_cells_per_perturbation=int(snk.params.max_cells_per_pert),
        random_seed=int(snk.params.seed),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
