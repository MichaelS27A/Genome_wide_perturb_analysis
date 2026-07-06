#!/usr/bin/env python3
"""PS downstream DE: select high-score cells per perturbation and run DE vs controls."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _to_numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        raise RuntimeError(f"Column '{col}' not found in input table.")
    return pd.to_numeric(df[col], errors="coerce")


def _select_cells(ps_df: pd.DataFrame, score_col: str, score_mode: str, score_quantile: float) -> pd.DataFrame:
    score_mode = str(score_mode).lower().strip()
    if score_mode not in {"top_positive", "top_negative", "top_absolute"}:
        raise RuntimeError("score_mode must be one of: top_positive, top_negative, top_absolute")

    out = []
    for pert, g in ps_df.groupby("perturbation", sort=False):
        s = pd.to_numeric(g[score_col], errors="coerce")
        valid = s.notna()
        if not bool(valid.any()):
            continue
        gv = g.loc[valid].copy()
        sv = s.loc[valid].to_numpy(dtype=float)

        if score_mode == "top_positive":
            metric = sv
        elif score_mode == "top_negative":
            metric = -sv
        else:
            metric = np.abs(sv)

        thr = float(np.quantile(metric, score_quantile))
        keep = metric >= thr
        if not bool(np.any(keep)):
            continue
        sel = gv.loc[keep].copy()
        sel["selection_metric"] = metric[keep]
        sel["selection_threshold"] = thr
        out.append(sel)

    if not out:
        return pd.DataFrame(columns=list(ps_df.columns) + ["selection_metric", "selection_threshold"])
    return pd.concat(out, axis=0, ignore_index=True)


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


def run_analysis(args: argparse.Namespace) -> None:
    args.h5ad = Path(args.h5ad)
    args.cell_scores = Path(args.cell_scores)
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if not (0.0 <= float(args.score_quantile) <= 1.0):
        raise RuntimeError("score_quantile must be in [0, 1].")

    ps = pd.read_csv(args.cell_scores, sep="\t", compression="infer", dtype="string")
    if ps.empty:
        raise RuntimeError("PS cell-score table is empty.")
    if "cell_barcode" not in ps.columns or "perturbation" not in ps.columns:
        raise RuntimeError("PS cell-score table must have columns: cell_barcode, perturbation.")

    ps["cell_barcode"] = ps["cell_barcode"].astype(str)
    ps["perturbation"] = ps["perturbation"].astype(str)
    ps[args.score_col] = _to_numeric_series(ps, args.score_col)
    ps = ps.dropna(subset=[args.score_col]).copy()
    ps = ps[ps["perturbation"] != args.control_label].copy()

    selected = _select_cells(
        ps_df=ps,
        score_col=args.score_col,
        score_mode=args.score_mode,
        score_quantile=float(args.score_quantile),
    )
    if selected.empty:
        raise RuntimeError("No cells selected from PS scores with current selection settings.")

    selected_counts = (
        selected.groupby("perturbation", as_index=False)
        .agg(n_selected_cells=("cell_barcode", "nunique"))
        .sort_values("n_selected_cells", ascending=False)
        .reset_index(drop=True)
    )
    keep_perts = set(
        selected_counts.loc[selected_counts["n_selected_cells"] >= int(args.min_selected_cells), "perturbation"]
        .astype(str)
        .tolist()
    )
    selected = selected[selected["perturbation"].isin(keep_perts)].copy()
    if selected.empty:
        raise RuntimeError(
            f"No perturbations with >= {args.min_selected_cells} selected cells after PS score filtering."
        )

    selected = selected.sort_values(["perturbation", args.score_col], ascending=[True, False]).drop_duplicates(
        subset=["cell_barcode"], keep="first"
    )
    selected_map = selected.set_index("cell_barcode")["perturbation"]

    adata_backed = sc.read_h5ad(args.h5ad, backed="r")
    obs = adata_backed.obs[[args.pert_col]].copy()
    obs.index = obs.index.astype(str)
    obs[args.pert_col] = obs[args.pert_col].astype(str)

    control_barcodes = obs.index[obs[args.pert_col] == args.control_label].to_numpy(dtype=str)
    if control_barcodes.size == 0:
        adata_backed.file.close()
        raise RuntimeError(f"No control cells found for label '{args.control_label}'.")

    if int(args.max_control_cells) > 0 and control_barcodes.size > int(args.max_control_cells):
        rng = np.random.default_rng(int(args.random_seed))
        control_barcodes = rng.choice(control_barcodes, size=int(args.max_control_cells), replace=False)

    keep_barcodes = pd.Index(selected_map.index.tolist() + control_barcodes.tolist()).unique()
    obs_names = pd.Index(adata_backed.obs_names.astype(str))
    keep_mask = obs_names.isin(keep_barcodes)
    if int(keep_mask.sum()) == 0:
        adata_backed.file.close()
        raise RuntimeError("No overlap between selected/control barcodes and adata.obs_names.")

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

    adata_de = adata.copy()
    sc.pp.normalize_total(adata_de)
    sc.pp.log1p(adata_de)

    de_counts = adata_de.obs["de_group"].value_counts()
    valid_groups = [
        g
        for g in de_counts.index.tolist()
        if g != args.control_label and int(de_counts[g]) >= int(args.min_selected_cells)
    ]
    if valid_groups:
        sc.tl.rank_genes_groups(
            adata_de,
            groupby="de_group",
            groups=valid_groups,
            reference=args.control_label,
            method="wilcoxon",
            pts=True,
            n_genes=min(int(args.n_top_de_genes), int(adata_de.n_vars)),
        )
        de_df = extract_rank_genes_groups(
            adata_de, groups=valid_groups, n_top=min(int(args.n_top_de_genes), int(adata_de.n_vars))
        )
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

    selected_out = selected[
        ["cell_barcode", "perturbation", args.score_col, "selection_metric", "selection_threshold"]
    ].rename(columns={args.score_col: "ps_score"})
    counts_out = (
        selected_out.groupby("perturbation", as_index=False)
        .agg(
            n_selected_cells=("cell_barcode", "nunique"),
            mean_ps_score=("ps_score", "mean"),
            median_ps_score=("ps_score", "median"),
        )
        .sort_values("n_selected_cells", ascending=False)
        .reset_index(drop=True)
    )

    selected_out.to_csv(args.outdir / "selected_ps_cells.tsv.gz", sep="\t", index=False, compression="gzip")
    counts_out.to_csv(args.outdir / "perturbation_selection_counts.tsv.gz", sep="\t", index=False, compression="gzip")
    de_df.to_csv(args.outdir / "perturbation_differential_genes.tsv.gz", sep="\t", index=False, compression="gzip")

    meta = {
        "h5ad": str(args.h5ad),
        "cell_scores": str(args.cell_scores),
        "pert_col": args.pert_col,
        "control_label": args.control_label,
        "score_col": args.score_col,
        "score_mode": args.score_mode,
        "score_quantile": float(args.score_quantile),
        "min_selected_cells": int(args.min_selected_cells),
        "max_control_cells": int(args.max_control_cells),
        "n_top_de_genes": int(args.n_top_de_genes),
        "random_seed": int(args.random_seed),
        "n_selected_cells": int(selected_out["cell_barcode"].nunique()),
        "n_selected_perturbations": int(counts_out.shape[0]),
        "n_de_rows": int(de_df.shape[0]),
    }
    (args.outdir / "de_meta.json").write_text(json.dumps(meta, indent=2))
    (args.outdir / "done.txt").write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--cell-scores", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--score-col", type=str, default="ps_score")
    ap.add_argument("--score-mode", type=str, default="top_positive")
    ap.add_argument("--score-quantile", type=float, default=0.90)
    ap.add_argument("--min-selected-cells", type=int, default=20)
    ap.add_argument("--max-control-cells", type=int, default=50000)
    ap.add_argument("--n-top-de-genes", type=int, default=200)
    ap.add_argument("--random-seed", type=int, default=0)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        h5ad=Path(str(snk.input.h5ad)),
        cell_scores=Path(str(snk.input.cell_scores)),
        outdir=Path(str(snk.params.outdir)),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        score_col=str(snk.params.score_col),
        score_mode=str(snk.params.score_mode),
        score_quantile=float(snk.params.score_quantile),
        min_selected_cells=int(snk.params.min_selected),
        max_control_cells=int(snk.params.max_control_cells),
        n_top_de_genes=int(snk.params.n_top_de_genes),
        random_seed=int(snk.params.random_seed),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
