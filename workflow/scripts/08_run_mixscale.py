#!/usr/bin/env python3
"""Run pertpy Mixscale method and write workflow-compatible outputs."""


import argparse
import json
import os
import platform
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

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import read_elem
from scipy import sparse


def _open_minimal_backed_adata(h5ad_path: Path) -> tuple[ad.AnnData, h5py.File | None]:
    """Open H5AD without eagerly materializing heavy optional groups (e.g. layers).

    Falls back to anndata backed reader for files whose X metadata does not expose
    the sparse encoding attributes expected by sparse_dataset().
    """
    h5 = h5py.File(h5ad_path, "r")
    x_obj = h5.get("X")
    has_sparse_attrs = (
        isinstance(x_obj, h5py.Group)
        and ("encoding-type" in x_obj.attrs or "h5sparse_format" in x_obj.attrs)
    )
    if not has_sparse_attrs:
        h5.close()
        return ad.read_h5ad(h5ad_path, backed="r"), None
    try:
        X = sparse_dataset(h5["X"])
        obs = read_elem(h5["obs"])
        var = read_elem(h5["var"])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        return adata, h5
    except Exception as e:
        h5.close()
        msg = str(e)
        if "encoding-type" not in msg and "h5sparse_format" not in msg:
            raise
        # Compatibility path for variant/legacy H5AD matrix encoding.
        adata = ad.read_h5ad(h5ad_path, backed="r")
        return adata, None


def _close_backing_resources(adata: ad.AnnData, h5_handle: h5py.File | None) -> None:
    if h5_handle is not None:
        try:
            h5_handle.close()
        except Exception:
            pass
        return
    file_handle = getattr(adata, "file", None)
    if file_handle is not None:
        try:
            file_handle.close()
        except Exception:
            pass


def _read_csc_nnz_per_gene(h5: h5py.File | None) -> tuple[np.ndarray | None, int]:
    """Read per-gene nnz for CSC X if available; returns (nnz_per_gene, bytes_per_nnz)."""
    if h5 is None:
        return None, 12
    try:
        xg = h5["X"]
        if "indptr" not in xg:
            return None, 12
        indptr = np.asarray(xg["indptr"][...], dtype=np.int64)
        if indptr.ndim != 1 or indptr.size < 2:
            return None, 12
        nnz = np.diff(indptr).astype(np.int64, copy=False)
        bytes_per_nnz = int(xg["data"].dtype.itemsize + xg["indices"].dtype.itemsize)
        if bytes_per_nnz <= 0:
            bytes_per_nnz = 12
        return nnz, bytes_per_nnz
    except Exception:
        return None, 12


def _select_csc_gene_subset(
    var: pd.DataFrame,
    nnz_per_gene: np.ndarray | None,
    max_genes: int,
    max_total_nnz: int,
) -> tuple[np.ndarray, dict[str, int]]:
    if max_genes <= 0 or max_genes >= int(var.shape[0]):
        full = np.arange(int(var.shape[0]), dtype=np.int64)
        return full, {
            "selected_genes": int(full.size),
            "selected_total_nnz": int(-1),
            "per_gene_nnz_cap": int(-1),
            "selector_used_budget": int(0),
        }

    n_vars = int(var.shape[0])
    if nnz_per_gene is None or int(nnz_per_gene.shape[0]) != n_vars:
        nnz_per_gene = np.zeros(n_vars, dtype=np.int64)
    nnz_per_gene = nnz_per_gene.astype(np.int64, copy=False)

    score_cols = (
        "highly_variable_rank",
        "dispersions_norm",
        "variance",
        "mean_counts",
        "total_counts",
        "n_cells_by_counts",
        "n_cells",
    )
    order = np.arange(n_vars, dtype=np.int64)
    for col in score_cols:
        if col not in var.columns:
            continue
        scores = pd.to_numeric(var[col], errors="coerce").to_numpy(dtype=np.float64, copy=False)
        valid = np.isfinite(scores)
        if not bool(valid.any()):
            continue
        if col == "highly_variable_rank":
            safe_scores = np.where(valid, scores, np.inf)
            order = np.argsort(safe_scores, kind="stable").astype(np.int64, copy=False)
        else:
            safe_scores = np.where(valid, scores, -np.inf)
            order = np.argsort(safe_scores, kind="stable")[::-1].astype(np.int64, copy=False)
        break

    nonzero_mask = nnz_per_gene > 0
    selected: list[int] = []
    selected_set: set[int] = set()
    selected_total_nnz = 0

    if max_total_nnz > 0:
        per_gene_cap = max(1, int(max_total_nnz) // max(1, int(max_genes)))
        min_candidates = max(100, int(max_genes) // 3)
        candidate_order = np.array([], dtype=np.int64)
        for fac in (1, 2, 4, 8, 16, 32, 64, 128):
            cap = int(per_gene_cap) * int(fac)
            cand_mask = nonzero_mask & (nnz_per_gene <= cap)
            candidate_order = order[cand_mask[order]]
            if int(candidate_order.size) >= min_candidates or fac == 128:
                per_gene_cap = cap
                break
    else:
        per_gene_cap = -1
        candidate_order = order[nonzero_mask[order]]

    for g in candidate_order:
        gi = int(g)
        gn = int(nnz_per_gene[gi])
        if max_total_nnz > 0 and selected_total_nnz + gn > int(max_total_nnz):
            continue
        selected.append(gi)
        selected_set.add(gi)
        selected_total_nnz += gn
        if len(selected) >= int(max_genes):
            break

    if len(selected) < int(max_genes):
        remainder = order[
            nonzero_mask[order]
            & np.array([int(x) not in selected_set for x in order], dtype=bool)
        ]
        for g in remainder:
            gi = int(g)
            gn = int(nnz_per_gene[gi])
            if max_total_nnz > 0 and selected_total_nnz + gn > int(max_total_nnz):
                continue
            selected.append(gi)
            selected_set.add(gi)
            selected_total_nnz += gn
            if len(selected) >= int(max_genes):
                break

    if not selected:
        fallback = np.arange(min(int(max_genes), n_vars), dtype=np.int64)
        return fallback, {
            "selected_genes": int(fallback.size),
            "selected_total_nnz": int(-1),
            "per_gene_nnz_cap": int(-1),
            "selector_used_budget": int(0),
        }

    out = np.sort(np.asarray(selected, dtype=np.int64))
    return out, {
        "selected_genes": int(out.size),
        "selected_total_nnz": int(selected_total_nnz),
        "per_gene_nnz_cap": int(per_gene_cap),
        "selector_used_budget": int(max_total_nnz > 0),
    }


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
    csc_max_genes: int = 1000,
    csc_max_total_nnz: int = 120000000,
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
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")
    adata.obs[args.pert_col] = adata.obs[args.pert_col].astype(str)
    if args.control_label not in set(adata.obs[args.pert_col].tolist()):
        raise RuntimeError(f"control label '{args.control_label}' not present in '{args.pert_col}'")

    sc.pp.normalize_total(adata, target_sum=args.normalize_target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(adata, use_highly_variable=bool(args.use_hvg_for_pca))

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
            "auto_batch_max_elements": int(args.auto_batch_max_elements),
            "auto_batch_size": int(args.auto_batch_size),
            "csc_max_genes": int(args.csc_max_genes),
            "csc_max_total_nnz": int(args.csc_max_total_nnz),
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
        default=1000,
        help=(
            "For CSC-backed H5AD matrices, prefilter to at most this many genes before row slicing "
            "to avoid large-memory backed slicing paths. Set 0 to disable."
        ),
    )
    ap.add_argument(
        "--csc-max-total-nnz",
        type=int,
        default=120000000,
        help=(
            "For CSC-backed H5AD matrices, cap total nonzeros selected during gene prefiltering. "
            "Set 0 to disable this nnz budget."
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
        mixscale_logfc_threshold=float(snk.params.logfc_threshold),
        mixscale_pval_cutoff=float(snk.params.pval_cutoff),
        mixscale_min_de_genes=int(snk.params.min_de_genes),
        mixscale_max_de_genes=int(snk.params.max_de_genes),
        batch_size=int(snk.params.batch_size),
        auto_batch_max_elements=int(snk.params.auto_batch_max_elements),
        auto_batch_size=int(snk.params.auto_batch_size),
        csc_max_genes=int(getattr(snk.params, "csc_max_genes", 1000)),
        csc_max_total_nnz=int(getattr(snk.params, "csc_max_total_nnz", 120000000)),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
