#!/usr/bin/env python3
"""Run pertpy Mixscale method and write workflow-compatible outputs."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def _load_by_barcodes(h5ad_path: Path, barcodes: list[str]) -> ad.AnnData:
    if not barcodes:
        raise RuntimeError("No chunk barcodes provided")

    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs_names = pd.Index(adata.obs_names.astype(str))
    keep = obs_names.isin(barcodes)
    if int(keep.sum()) == 0:
        adata.file.close()
        raise RuntimeError("No overlap between chunk barcodes and adata.obs_names")

    try:
        try:
            return adata[keep].to_memory()
        except Exception:
            idx = np.where(keep)[0]
            parts: list[ad.AnnData] = []
            batch = 512
            for i in range(0, len(idx), batch):
                j = min(i + batch, len(idx))
                parts.append(adata[idx[i:j]].to_memory())
            if not parts:
                raise RuntimeError("Chunk subset fallback failed: no cells loaded")
            return ad.concat(parts, axis=0, join="outer", merge="same")
    finally:
        adata.file.close()


def _load_full_with_optional_subsample(h5ad_path: Path, max_cells: int, seed: int) -> ad.AnnData:
    adata = sc.read_h5ad(h5ad_path, backed="r")
    n_obs = int(adata.n_obs)
    take_idx: np.ndarray | None = None
    if max_cells > 0 and n_obs > max_cells:
        rng = np.random.default_rng(seed)
        take_idx = np.sort(rng.choice(n_obs, size=max_cells, replace=False))
    try:
        if take_idx is None:
            return adata.to_memory()
        try:
            return adata[take_idx].to_memory()
        except Exception:
            parts: list[ad.AnnData] = []
            batch = 512
            for i in range(0, len(take_idx), batch):
                j = min(i + batch, len(take_idx))
                parts.append(adata[take_idx[i:j]].to_memory())
            if not parts:
                raise RuntimeError("Subsampling failed: no cells loaded")
            return ad.concat(parts, axis=0, join="outer", merge="same")
    finally:
        adata.file.close()


def load_h5ad(h5ad_path: Path, max_cells: int, seed: int, chunk_cells: Path | None = None) -> ad.AnnData:
    if chunk_cells is not None:
        df = pd.read_csv(chunk_cells, sep="\t", compression="infer", dtype="string")
        barcodes = df["cell_barcode"].dropna().astype(str).unique().tolist()
        return _load_by_barcodes(h5ad_path=h5ad_path, barcodes=barcodes)
    return _load_full_with_optional_subsample(h5ad_path=h5ad_path, max_cells=max_cells, seed=seed)


def choose_targets(obs: pd.DataFrame, pert_col: str, control_label: str, min_cells: int, max_perturbations: int) -> list[str]:
    counts = obs[pert_col].astype(str).value_counts()
    targets = [p for p, c in counts.items() if p != control_label and int(c) >= int(min_cells)]
    if max_perturbations > 0:
        targets = targets[:max_perturbations]
    return targets


def build_de_table(adata: ad.AnnData, targets: list[str], max_de_genes: int) -> pd.DataFrame:
    de_map = adata.uns.get("mixscale_de_genes", {})
    rows: list[dict[str, object]] = []
    for pert in targets:
        genes = de_map.get(pert, [])
        genes = np.asarray(genes).astype(str).tolist()
        if max_de_genes > 0:
            genes = genes[:max_de_genes]
        for rank, gene in enumerate(genes, start=1):
            rows.append(
                {
                    "perturbation": pert,
                    "gene": gene,
                    "rank": rank,
                    "p_weight": np.nan,
                }
            )
    cols = ["perturbation", "gene", "rank", "p_weight"]
    return pd.DataFrame(rows, columns=cols)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--chunk-cells", type=Path, default=None)
    ap.add_argument("--chunk-id", type=str, default=None)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--min-cells-per-perturbation", type=int, default=30)
    ap.add_argument("--max-perturbations", type=int, default=0)
    ap.add_argument("--max-cells", type=int, default=0)
    ap.add_argument("--random-seed", type=int, default=0)
    ap.add_argument("--use-hvg-for-pca", action="store_true")
    ap.add_argument("--normalize-target-sum", type=float, default=1e4)
    ap.add_argument("--mixscale-logfc-threshold", type=float, default=0.10)
    ap.add_argument("--mixscale-pval-cutoff", type=float, default=0.05)
    ap.add_argument("--mixscale-min-de-genes", type=int, default=5)
    ap.add_argument("--mixscale-max-de-genes", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=0)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = load_h5ad(
        args.h5ad,
        max_cells=args.max_cells,
        seed=args.random_seed,
        chunk_cells=args.chunk_cells,
    )
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")
    adata.obs[args.pert_col] = adata.obs[args.pert_col].astype(str)
    if args.control_label not in set(adata.obs[args.pert_col].tolist()):
        raise RuntimeError(f"control label '{args.control_label}' not present in '{args.pert_col}'")

    sc.pp.normalize_total(adata, target_sum=args.normalize_target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(adata, use_highly_variable=bool(args.use_hvg_for_pca))

    import pertpy as pt

    if not hasattr(pt.tl, "Mixscale"):
        raise RuntimeError(
            "This pertpy build does not provide Mixscale. "
            "Update to a newer pertpy release (for example >=1.1.1)."
        )

    ms = pt.tl.Mixscale()
    batch_size = args.batch_size if args.batch_size and args.batch_size > 0 else None
    ms.perturbation_signature(
        adata=adata,
        pert_key=args.pert_col,
        control=args.control_label,
        split_by=None,
        batch_size=batch_size,
    )
    ms.mixscale(
        adata=adata,
        pert_key=args.pert_col,
        control=args.control_label,
        new_class_name="mixscale_score",
        layer="X_pert",
        min_de_genes=args.mixscale_min_de_genes,
        max_de_genes=args.mixscale_max_de_genes,
        logfc_threshold=args.mixscale_logfc_threshold,
        pval_cutoff=args.mixscale_pval_cutoff,
        split_by=None,
        random_state=args.random_seed,
    )

    if "mixscale_score" not in adata.obs.columns:
        raise RuntimeError("mixscale_score not found in adata.obs after Mixscale run")

    cell_scores = pd.DataFrame(
        {
            "cell_barcode": adata.obs_names.astype(str),
            "perturbation": adata.obs[args.pert_col].astype(str).values,
            "mixscale_score": pd.to_numeric(adata.obs["mixscale_score"], errors="coerce").values,
        }
    )

    targets = choose_targets(
        adata.obs,
        pert_col=args.pert_col,
        control_label=args.control_label,
        min_cells=args.min_cells_per_perturbation,
        max_perturbations=args.max_perturbations,
    )
    de_tbl = build_de_table(adata, targets=targets, max_de_genes=args.mixscale_max_de_genes)

    cell_scores.to_csv(args.outdir / "cell_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    de_tbl.to_csv(args.outdir / "perturbation_de.tsv.gz", sep="\t", index=False, compression="gzip")

    meta = {
        "method": "Mixscale_pertpy",
        "h5ad": str(args.h5ad),
        "pert_col": args.pert_col,
        "control_label": args.control_label,
        "n_cells": int(adata.n_obs),
        "n_perturbations_total": int(adata.obs[args.pert_col].nunique()),
        "n_perturbations_tested": int(len(targets)),
        "max_cells": int(args.max_cells),
        "chunk_cells": str(args.chunk_cells) if args.chunk_cells else None,
        "chunk_id": args.chunk_id,
        "parameters": {
            "use_hvg_for_pca": bool(args.use_hvg_for_pca),
            "normalize_target_sum": float(args.normalize_target_sum),
            "mixscale_logfc_threshold": float(args.mixscale_logfc_threshold),
            "mixscale_pval_cutoff": float(args.mixscale_pval_cutoff),
            "mixscale_min_de_genes": int(args.mixscale_min_de_genes),
            "mixscale_max_de_genes": int(args.mixscale_max_de_genes),
            "batch_size": batch_size,
            "min_cells_per_perturbation": int(args.min_cells_per_perturbation),
            "max_perturbations": int(args.max_perturbations),
            "random_seed": int(args.random_seed),
        },
        "package_versions": {
            "python": platform.python_version(),
            "anndata": ad.__version__,
            "scanpy": sc.__version__,
            "pertpy": pt.__version__,
        },
    }
    (args.outdir / "method_meta.json").write_text(json.dumps(meta, indent=2))
    (args.outdir / "done.txt").write_text("ok\n")


if __name__ == "__main__":
    main()
