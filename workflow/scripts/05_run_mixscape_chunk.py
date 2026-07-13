#!/usr/bin/env python3
"""Run chunk-level Mixscape and produce perturbation effect summaries."""


import argparse
import inspect
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

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import read_elem
from scipy import stats
from scipy import sparse


def detect_perturbed_mask(obs: pd.DataFrame) -> pd.Series:
    """Heuristic for perturbed-class cells across different Mixscape label variants."""
    candidates = [c for c in obs.columns if "mixscape" in c.lower() and "class" in c.lower()]
    if not candidates:
        return pd.Series(False, index=obs.index)

    preferred = [c for c in ("mixscape_class_global", "mixscape_class") if c in candidates]
    remaining = [c for c in candidates if c not in preferred]
    ordered = preferred + remaining

    for col in ordered:
        vals = obs[col].astype(str).str.lower()
        mask = vals.str.contains("ko|perturbed", regex=True, na=False)
        if bool(mask.any()):
            return mask
    return pd.Series(False, index=obs.index)


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


def load_subset(
    h5ad_path: Path,
    chunk_cells_tsv: Path,
    csc_max_genes: int,
    csc_max_total_nnz: int,
) -> ad.AnnData:
    chunk_df = pd.read_csv(chunk_cells_tsv, sep="\t", compression="infer", dtype="string")
    barcodes = chunk_df["cell_barcode"].dropna().astype(str).unique().tolist()

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
            "[mixscape] CSC-backed matrix detected; prefiltering genes before row slicing "
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
        sub = work[keep].to_memory()
    except Exception as e:
        idx = np.where(keep)[0]
        batch = 512
        parts: list[ad.AnnData] = []
        for i in range(0, len(idx), batch):
            j = min(i + batch, len(idx))
            parts.append(work[idx[i:j]].to_memory())
        if not parts:
            raise RuntimeError("Subset fallback failed: no cells loaded") from e
        sub = ad.concat(parts, axis=0, join="outer", merge="same")
    finally:
        _close_backing_resources(adata, h5_handle)

    chunk_meta = chunk_df.drop_duplicates(subset=["cell_barcode"]).set_index("cell_barcode")
    common = sub.obs_names.intersection(chunk_meta.index)
    sub.obs = sub.obs.copy()
    sub.obs.loc[common, "chunk_gene_target"] = chunk_meta.loc[common, "gene_target"].astype(str).values
    sub.obs.loc[common, "chunk_is_control"] = chunk_meta.loc[common, "is_control"].astype(int).values
    return sub


def run_mixscape(
    adata: ad.AnnData,
    pert_col: str,
    control_label: str,
    use_hvg_for_pca: bool = False,
    batch_size: int | None = None,
    auto_batch_max_elements: int = 800_000_000,
    auto_batch_size: int = 2000,
    normalize_target_sum: float = 1e4,
    mixscape_logfc_threshold: float = 0.10,
    mixscape_pval_cutoff: float = 0.05,
) -> None:
    import pertpy as pt

    sc.pp.normalize_total(adata, target_sum=normalize_target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(adata, use_highly_variable=use_hvg_for_pca)

    # Sparse perturbation_signature code paths in pertpy have shown instability
    # on very wide matrices; operate on dense float32 for robustness.
    if sparse.issparse(adata.X):
        print(
            f"[mixscape] Converting sparse matrix to dense float32 "
            f"(cells={adata.n_obs}, genes={adata.n_vars})",
            flush=True,
        )
        adata.X = adata.X.toarray().astype(np.float32, copy=False)

    matrix_elements = int(adata.n_obs) * int(adata.n_vars)
    if (
        batch_size is None
        and auto_batch_max_elements > 0
        and matrix_elements > int(auto_batch_max_elements)
    ):
        batch_size = max(1, int(auto_batch_size))
        print(
            "[mixscape] Auto-enabling batching due to large matrix: "
            f"elements={matrix_elements} threshold={int(auto_batch_max_elements)} "
            f"batch_size={batch_size}",
            flush=True,
        )

    ms = pt.tl.Mixscape()
    split_by = None
    ms.perturbation_signature(
        adata,
        pert_key=pert_col,
        control=control_label,
        split_by=split_by,
        batch_size=batch_size,
    )
    mixscape_kwargs = dict(
        adata=adata,
        control=control_label,
        layer="X_pert",
        logfc_threshold=mixscape_logfc_threshold,
        pval_cutoff=mixscape_pval_cutoff,
    )
    if "pert_key" in inspect.signature(ms.mixscape).parameters:
        mixscape_kwargs["pert_key"] = pert_col
    else:
        mixscape_kwargs["labels"] = pert_col
    ms.mixscape(**mixscape_kwargs)


def summarize_chunk(
    adata: ad.AnnData,
    pert_col: str,
    control_label: str,
    pca_dims: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    obs = adata.obs.copy()
    obs[pert_col] = obs[pert_col].astype(str)

    if "X_pca" not in adata.obsm:
        raise RuntimeError("X_pca not found; expected PCA run before summarization")

    pca = np.asarray(adata.obsm["X_pca"])
    nd = min(pca_dims, pca.shape[1])
    pca = pca[:, :nd]

    ctrl_mask = obs[pert_col] == control_label
    if int(ctrl_mask.sum()) == 0:
        raise RuntimeError(f"No control cells found for label '{control_label}' in this chunk")

    ctrl_vec = pca[ctrl_mask.values].mean(axis=0)
    perturbed_mask = detect_perturbed_mask(obs)

    stats_rows = []
    effect_rows = []
    labels_rows = []

    out_cols = [pert_col]
    mix_cols = [c for c in obs.columns if "mixscape" in c.lower()]
    out_cols.extend(mix_cols)
    lbl = obs[out_cols].copy()
    lbl.insert(0, "cell_barcode", lbl.index.astype(str))
    lbl["is_predicted_perturbed"] = perturbed_mask.values.astype(int)
    labels_rows.append(lbl)

    for pert in sorted(obs[pert_col].unique().tolist()):
        mask = obs[pert_col] == pert
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue

        vec = pca[mask.values]
        mean_vec = vec.mean(axis=0)
        delta = mean_vec - ctrl_vec
        l2 = float(np.linalg.norm(delta))

        if pert == control_label:
            t_stat = np.nan
            p_val = np.nan
            n_pred = int(perturbed_mask[mask].sum())
        else:
            t_res = stats.ttest_ind(vec[:, 0], pca[ctrl_mask.values][:, 0], equal_var=False, nan_policy="omit")
            t_stat = float(t_res.statistic) if np.isfinite(t_res.statistic) else np.nan
            p_val = float(t_res.pvalue) if np.isfinite(t_res.pvalue) else np.nan
            n_pred = int(perturbed_mask[mask].sum())

        stats_rows.append(
            {
                "perturbation": pert,
                "n_cells": n_cells,
                "n_predicted_perturbed": n_pred,
                "frac_predicted_perturbed": (n_pred / n_cells) if n_cells else np.nan,
                "pc1_t_stat_vs_control": t_stat,
                "pc1_p_value_vs_control": p_val,
                "effect_l2_vs_control": l2,
            }
        )

        effect_rec = {"perturbation": pert, "n_cells": n_cells, "effect_l2_vs_control": l2}
        for i, x in enumerate(delta, start=1):
            effect_rec[f"delta_pc{i}"] = float(x)
        effect_rows.append(effect_rec)

    stats_df = pd.DataFrame(stats_rows).sort_values(["perturbation"]).reset_index(drop=True)
    effects_df = pd.DataFrame(effect_rows).sort_values(["perturbation"]).reset_index(drop=True)
    labels_df = pd.concat(labels_rows, axis=0)

    return stats_df, effects_df, labels_df


def run_analysis(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_subset(
        args.h5ad,
        args.chunk_cells,
        csc_max_genes=args.csc_max_genes,
        csc_max_total_nnz=args.csc_max_total_nnz,
    )
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")

    batch_size = args.batch_size if args.batch_size and args.batch_size > 0 else None

    run_mixscape(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        use_hvg_for_pca=args.use_hvg_for_pca,
        batch_size=batch_size,
        auto_batch_max_elements=args.auto_batch_max_elements,
        auto_batch_size=args.auto_batch_size,
        normalize_target_sum=args.normalize_target_sum,
        mixscape_logfc_threshold=args.mixscape_logfc_threshold,
        mixscape_pval_cutoff=args.mixscape_pval_cutoff,
    )
    stats_df, effects_df, labels_df = summarize_chunk(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        pca_dims=args.pca_dims,
    )

    sel = labels_df.copy()
    sel["perturbation"] = sel[args.pert_col].astype(str)
    sel = sel[
        (sel["is_predicted_perturbed"].astype(int) == 1) & (sel["perturbation"] != args.control_label)
    ][["cell_barcode", "perturbation"]].drop_duplicates()
    sel["chunk_id"] = str(args.chunk_id) if args.chunk_id is not None else args.chunk_cells.stem

    stats_df.to_csv(args.output_dir / "perturbation_stats.tsv", sep="\t", index=False)
    effects_df.to_csv(args.output_dir / "perturbation_effects_pca.tsv", sep="\t", index=False)
    labels_df.to_csv(args.output_dir / "cell_mixscape_labels.tsv.gz", sep="\t", index=False, compression="gzip")
    sel.to_csv(args.output_dir / "selected_perturbed_cells.tsv.gz", sep="\t", index=False, compression="gzip")

    if args.write_subset:
        adata.write_h5ad(args.output_dir / "subset_annotated.h5ad", compression="gzip")

    run_meta = {
        "h5ad": str(args.h5ad),
        "chunk_cells": str(args.chunk_cells),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "pert_col": args.pert_col,
        "control_label": args.control_label,
        "pca_dims": args.pca_dims,
        "use_hvg_for_pca": bool(args.use_hvg_for_pca),
        "batch_size": batch_size,
        "auto_batch_max_elements": int(args.auto_batch_max_elements),
        "auto_batch_size": int(args.auto_batch_size),
        "csc_max_genes": int(args.csc_max_genes),
        "csc_max_total_nnz": int(args.csc_max_total_nnz),
        "normalize_target_sum": args.normalize_target_sum,
        "mixscape_logfc_threshold": args.mixscape_logfc_threshold,
        "mixscape_pval_cutoff": args.mixscape_pval_cutoff,
        "n_selected_perturbed_cells": int(sel.shape[0]),
    }
    (args.output_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
    (args.output_dir / "done.txt").write_text("ok\n")

    print(f"Wrote {args.output_dir / 'done.txt'}")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--chunk-cells", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--pca-dims", type=int, default=20)
    ap.add_argument("--chunk-id", type=str, default=None)
    ap.add_argument(
        "--batch-size",
        type=int,
        default=0,
        help="Batch size for Mixscape perturbation_signature. Use 0 for full-chunk mode (no internal batching).",
    )
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
    ap.add_argument("--write-subset", action="store_true")
    ap.add_argument(
        "--use-hvg-for-pca",
        action="store_true",
        help="If set, PCA is computed on highly variable genes only. Default uses full transcriptome.",
    )
    ap.add_argument("--normalize-target-sum", type=float, default=1e4)
    ap.add_argument("--mixscape-logfc-threshold", type=float, default=0.10)
    ap.add_argument("--mixscape-pval-cutoff", type=float, default=0.05)
    return ap.parse_args()

def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        h5ad=Path(str(snk.input.h5ad)),
        chunk_cells=Path(str(snk.input.chunk_cells)),
        output_dir=Path(str(snk.params.outdir)),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        pca_dims=int(snk.params.pca_dims),
        chunk_id=str(getattr(snk.wildcards, "chunk", "")),
        batch_size=int(snk.params.batch_size),
        auto_batch_max_elements=int(snk.params.auto_batch_max_elements),
        auto_batch_size=int(snk.params.auto_batch_size),
        csc_max_genes=int(getattr(snk.params, "csc_max_genes", 0)),
        csc_max_total_nnz=int(getattr(snk.params, "csc_max_total_nnz", 0)),
        write_subset=bool(snk.params.write_subset),
        use_hvg_for_pca=bool(snk.params.use_hvg_for_pca),
        normalize_target_sum=float(snk.params.normalize_target_sum),
        mixscape_logfc_threshold=float(snk.params.logfc_threshold),
        mixscape_pval_cutoff=float(snk.params.pval_cutoff),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
