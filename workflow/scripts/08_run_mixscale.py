#!/usr/bin/env python3
"""Run pertpy Mixscale method and write workflow-compatible outputs."""


import argparse
import json
import os
import platform
from pathlib import Path

import sys
sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).resolve().parent))
from workflow_support import (
    _open_minimal_backed_adata, _close_backing_resources, _read_csc_nnz_per_gene,
    _select_csc_gene_subset, require_unique, require_barcodes, runtime_metadata, read_tables, require_raw_counts,
)


# Apply conservative defaults before importing numerical libraries.
for _env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_env, "1")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _load_by_barcodes(
    h5ad_path: Path,
    barcodes: list[str],
    csc_max_genes: int,
    csc_max_total_nnz: int,
) -> ad.AnnData:
    if not barcodes:
        raise RuntimeError("No chunk barcodes provided")

    adata, h5_handle = _open_minimal_backed_adata(h5ad_path)
    obs_names = pd.Index(adata.obs_names.astype(str))
    try:
        require_barcodes(obs_names, barcodes)
    except Exception:
        _close_backing_resources(adata, h5_handle)
        raise
    keep = obs_names.isin(barcodes)
    if int(keep.sum()) == 0:
        _close_backing_resources(adata, h5_handle)
        raise RuntimeError("No overlap between chunk barcodes and adata.obs_names")

    work = adata
    x_format = str(getattr(adata.X, "format", "")).lower()
    if x_format == "csc" and int(csc_max_genes) > 0 and int(csc_max_genes) < int(adata.n_vars):
        nnz_per_gene, bytes_per_nnz = _read_csc_nnz_per_gene(h5_handle)
        gene_idx, sel_meta = _select_csc_gene_subset(
            adata.var,
            nnz_per_gene=nnz_per_gene,
            max_genes=int(csc_max_genes),
            max_total_nnz=int(csc_max_total_nnz),
        )
        selected_total_nnz = int(sel_meta["selected_total_nnz"])
        est_gib = (
            (selected_total_nnz * bytes_per_nnz) / (1024.0**3)
            if selected_total_nnz >= 0
            else float("nan")
        )
        print(
            "[mixscale] CSC-backed matrix detected; prefiltering genes before row slicing "
            f"({len(gene_idx)} of {adata.n_vars}). "
            f"selection_total_nnz={selected_total_nnz} "
            f"per_gene_nnz_cap={int(sel_meta['per_gene_nnz_cap'])} "
            f"est_sparse_payload_gib={est_gib:.2f}",
            flush=True,
        )
        # Materialize the reduced-gene matrix first. Row slicing on backed CSC data
        # can otherwise trigger full-matrix reads and OOM.
        work = adata[:, gene_idx].to_memory()

    try:
        try:
            return work[keep].to_memory()
        except Exception:
            idx = np.where(keep)[0]
            parts: list[ad.AnnData] = []
            batch = 512
            for i in range(0, len(idx), batch):
                j = min(i + batch, len(idx))
                parts.append(work[idx[i:j]].to_memory())
            if not parts:
                raise RuntimeError("Chunk subset fallback failed: no cells loaded")
            return ad.concat(parts, axis=0, join="outer", merge="same")
    finally:
        _close_backing_resources(adata, h5_handle)


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


def load_h5ad(
    h5ad_path: Path,
    max_cells: int,
    seed: int,
    chunk_cells: Path | None = None,
    csc_max_genes: int = 0,
    csc_max_total_nnz: int = 0,
) -> ad.AnnData:
    if chunk_cells is not None:
        df = pd.read_csv(chunk_cells, sep="\t", compression="infer", dtype="string")
        barcodes = df["cell_barcode"].dropna().astype(str).unique().tolist()
        return _load_by_barcodes(
            h5ad_path=h5ad_path,
            barcodes=barcodes,
            csc_max_genes=csc_max_genes,
            csc_max_total_nnz=csc_max_total_nnz,
        )
    return _load_full_with_optional_subsample(h5ad_path=h5ad_path, max_cells=max_cells, seed=seed)


def attach_chunk_labels(adata: ad.AnnData, chunk_cells: Path, pert_col: str) -> None:
    """Attach the guide-call labels that selected this chunk to expression cells."""
    table = pd.read_csv(chunk_cells, sep="\t", compression="infer", dtype="string")
    required = {"cell_barcode", "gene_target"}
    if not required.issubset(table.columns):
        raise RuntimeError(f"Chunk table lacks required columns: {sorted(required)}")
    conflicts = table.groupby("cell_barcode", dropna=False)["gene_target"].nunique()
    if (conflicts > 1).any():
        raise RuntimeError("Chunk table assigns multiple perturbation labels to one cell")
    require_unique(table, ["cell_barcode"], "chunk cells")
    labels = table.set_index("cell_barcode")["gene_target"]
    labels.index = labels.index.astype(str)
    aligned = labels.reindex(adata.obs_names.astype(str))
    if aligned.isna().any():
        raise RuntimeError("Missing chunk-derived perturbation labels for expression cells")
    adata.obs = adata.obs.copy()
    adata.obs[pert_col] = aligned.astype(str).to_numpy()


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


def run_analysis(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)

    adata = load_h5ad(
        args.h5ad,
        max_cells=args.max_cells,
        seed=args.random_seed,
        chunk_cells=args.chunk_cells,
        csc_max_genes=args.csc_max_genes,
        csc_max_total_nnz=args.csc_max_total_nnz,
    )
    if args.chunk_cells is not None:
        attach_chunk_labels(adata, args.chunk_cells, args.pert_col)
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")
    adata.obs[args.pert_col] = adata.obs[args.pert_col].astype(str)
    if args.control_label not in set(adata.obs[args.pert_col].tolist()):
        raise RuntimeError(f"control label '{args.control_label}' not present in '{args.pert_col}'")

    require_raw_counts(adata.X)
    sc.pp.normalize_total(adata, target_sum=args.normalize_target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(
        adata,
        n_comps=min(int(args.pca_n_comps), max(1, min(adata.n_obs, adata.n_vars) - 1)),
        mask_var="highly_variable" if args.use_hvg_for_pca else None,
        random_state=args.random_seed,
    )

    # Sparse perturbation_signature code paths in pertpy are fragile for wide
    # matrices; convert to dense float32 before calling pertpy methods.
    if sparse.issparse(adata.X):
        print(
            f"[mixscale] Converting sparse matrix to dense float32 "
            f"(cells={adata.n_obs}, genes={adata.n_vars})",
            flush=True,
        )
        adata.X = adata.X.toarray().astype(np.float32, copy=False)

    import pertpy as pt

    if not hasattr(pt.tl, "Mixscale"):
        raise RuntimeError(
            "This pertpy build does not provide Mixscale. "
            "Update to a newer pertpy release (for example >=1.1.1)."
        )

    ms = pt.tl.Mixscale()
    batch_size = args.batch_size if args.batch_size and args.batch_size > 0 else None
    matrix_elements = int(adata.n_obs) * int(adata.n_vars)
    if (
        batch_size is None
        and args.auto_batch_max_elements > 0
        and matrix_elements > int(args.auto_batch_max_elements)
    ):
        batch_size = max(1, int(args.auto_batch_size))
        print(
            "[mixscale] Auto-enabling batching due to large matrix: "
            f"elements={matrix_elements} threshold={int(args.auto_batch_max_elements)} "
            f"batch_size={batch_size}",
            flush=True,
        )
    ms.perturbation_signature(
        adata=adata,
        pert_key=args.pert_col,
        control=args.control_label,
        ref_selection_mode=args.signature_ref_selection_mode,
        split_by=args.signature_split_by or None,
        n_neighbors=args.signature_n_neighbors,
        use_rep="X_pca",
        n_dims=min(args.signature_n_dims, int(adata.obsm["X_pca"].shape[1])),
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
        **runtime_metadata(adata),
        "effective_pca_n_comps": int(adata.obsm["X_pca"].shape[1]),
        "effective_signature_n_dims": min(args.signature_n_dims, int(adata.obsm["X_pca"].shape[1])),
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
            "pca_n_comps": int(args.pca_n_comps),
            "signature_n_neighbors": int(args.signature_n_neighbors),
            "signature_n_dims": int(args.signature_n_dims),
            "signature_ref_selection_mode": args.signature_ref_selection_mode,
            "signature_split_by": args.signature_split_by or None,
            "normalize_target_sum": float(args.normalize_target_sum),
            "mixscale_logfc_threshold": float(args.mixscale_logfc_threshold),
            "mixscale_pval_cutoff": float(args.mixscale_pval_cutoff),
            "mixscale_min_de_genes": int(args.mixscale_min_de_genes),
            "mixscale_max_de_genes": int(args.mixscale_max_de_genes),
            "batch_size": batch_size,
            "auto_batch_max_elements": int(args.auto_batch_max_elements),
            "auto_batch_size": int(args.auto_batch_size),
            "csc_max_genes": int(args.csc_max_genes),
            "csc_max_total_nnz": int(args.csc_max_total_nnz),
            "min_cells_per_perturbation": int(args.min_cells_per_perturbation),
            "max_perturbations": int(args.max_perturbations),
            "random_seed": int(args.random_seed),
        },
    }
    (args.outdir / "method_meta.json").write_text(json.dumps(meta, indent=2))
    (args.outdir / "done.txt").write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
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
    ap.add_argument("--pca-n-comps", type=int, default=60)
    ap.add_argument("--signature-n-neighbors", type=int, default=30)
    ap.add_argument("--signature-n-dims", type=int, default=60)
    ap.add_argument("--signature-ref-selection-mode", choices=("nn", "split_by"), default="nn")
    ap.add_argument("--signature-split-by", default="")
    ap.add_argument("--mixscale-logfc-threshold", type=float, default=0.10)
    ap.add_argument("--mixscale-pval-cutoff", type=float, default=0.05)
    ap.add_argument("--mixscale-min-de-genes", type=int, default=5)
    ap.add_argument("--mixscale-max-de-genes", type=int, default=100)
    ap.add_argument("--batch-size", type=int, default=0)
    ap.add_argument(
        "--auto-batch-max-elements",
        type=int,
        default=800_000_000,
        help=(
            "Auto-enable internal batching when n_cells*n_genes exceeds this threshold. "
            "Set 0 to disable auto-batching."
        ),
    )
    ap.add_argument(
        "--auto-batch-size",
        type=int,
        default=2000,
        help="Batch size to use when auto-batching is triggered.",
    )
    ap.add_argument(
        "--csc-max-genes",
        type=int,
        default=0,
        help=(
            "For CSC-backed H5AD matrices, prefilter to at most this many genes before row slicing "
            "to avoid large-memory backed slicing paths. Set >0 to enable."
        ),
    )
    ap.add_argument(
        "--csc-max-total-nnz",
        type=int,
        default=0,
        help=(
            "For CSC-backed H5AD matrices, cap total nonzeros selected during gene prefiltering. "
            "Set >0 to enable this nnz budget."
        ),
    )
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        h5ad=Path(str(snk.input.h5ad)),
        outdir=Path(str(snk.params.outdir)),
        chunk_cells=Path(str(snk.input.chunk_cells)),
        chunk_id=str(getattr(snk.wildcards, "chunk", "")),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        min_cells_per_perturbation=int(snk.params.min_cells),
        max_perturbations=int(snk.params.max_pert),
        max_cells=int(snk.params.max_cells),
        random_seed=int(snk.params.seed),
        use_hvg_for_pca=bool(snk.params.use_hvg_for_pca),
        normalize_target_sum=float(snk.params.normalize_target_sum),
        pca_n_comps=int(snk.params.pca_n_comps),
        signature_n_neighbors=int(snk.params.signature_n_neighbors),
        signature_n_dims=int(snk.params.signature_n_dims),
        signature_ref_selection_mode=str(snk.params.signature_ref_selection_mode),
        signature_split_by=str(snk.params.signature_split_by),
        mixscale_logfc_threshold=float(snk.params.logfc_threshold),
        mixscale_pval_cutoff=float(snk.params.pval_cutoff),
        mixscale_min_de_genes=int(snk.params.min_de_genes),
        mixscale_max_de_genes=int(snk.params.max_de_genes),
        batch_size=int(snk.params.batch_size),
        auto_batch_max_elements=int(snk.params.auto_batch_max_elements),
        auto_batch_size=int(snk.params.auto_batch_size),
        csc_max_genes=int(getattr(snk.params, "csc_max_genes", 0)),
        csc_max_total_nnz=int(getattr(snk.params, "csc_max_total_nnz", 0)),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
