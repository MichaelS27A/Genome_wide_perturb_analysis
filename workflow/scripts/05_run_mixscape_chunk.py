#!/usr/bin/env python3
"""Run chunk-level Mixscape and produce perturbation effect summaries."""


import argparse
import inspect
import json
import os
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
from scipy import stats
from scipy import sparse


def detect_perturbed_mask(obs: pd.DataFrame, perturbed_class: str = "KO") -> pd.Series:
    """Use pertpy 1.1.1's global classification contract."""
    if "mixscape_class_global" not in obs:
        raise RuntimeError("pertpy did not produce mixscape_class_global")
    return obs["mixscape_class_global"].astype(str).eq(perturbed_class)


def load_subset(
    h5ad_path: Path,
    chunk_cells_tsv: Path,
    pert_col: str,
    csc_max_genes: int,
    csc_max_total_nnz: int,
) -> ad.AnnData:
    chunk_df = pd.read_csv(chunk_cells_tsv, sep="\t", compression="infer", dtype="string")
    require_unique(chunk_df, ["cell_barcode"], "chunk cells")
    if chunk_df["gene_target"].isna().any():
        raise ValueError("Chunk contains missing perturbation labels")
    barcodes = chunk_df["cell_barcode"].astype(str).tolist()

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
    labels = chunk_meta["gene_target"].reindex(sub.obs_names.astype(str))
    if labels.isna().any():
        raise RuntimeError(f"Missing chunk-derived perturbation labels in '{pert_col}'")
    sub.obs[pert_col] = labels.astype(str).to_numpy()
    return sub


def run_mixscape(
    adata: ad.AnnData,
    pert_col: str,
    control_label: str,
    use_hvg_for_pca: bool = False,
    pca_n_comps: int = 60,
    signature_n_neighbors: int = 30,
    signature_n_dims: int = 60,
    signature_ref_selection_mode: str = "nn",
    signature_split_by: str | None = None,
    batch_size: int | None = None,
    auto_batch_max_elements: int = 800_000_000,
    auto_batch_size: int = 2000,
    normalize_target_sum: float = 1e4,
    mixscape_logfc_threshold: float = 0.10,
    mixscape_pval_cutoff: float = 0.05,
    mixscape_min_de_genes: int = 1,
    mixscape_test_method: str = "wilcoxon",
    mixscape_iter_num: int = 10,
    mixscape_scale: bool = True,
    mixscape_split_by: str | None = None,
    mixscape_perturbation_type: str = "KO",
    mixscape_random_state: int = 0,
) -> dict:
    import pertpy as pt

    require_raw_counts(adata.X)
    sc.pp.normalize_total(adata, target_sum=normalize_target_sum)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, subset=False)
    n_pca_genes = int(adata.var["highly_variable"].sum()) if use_hvg_for_pca else adata.n_vars
    max_pcs = min(int(pca_n_comps), adata.n_obs - 1, n_pca_genes - 1)
    if max_pcs < 1:
        raise ValueError("PCA needs at least two cells and two selected genes")
    sc.pp.pca(adata, n_comps=max_pcs, mask_var="highly_variable" if use_hvg_for_pca else None, random_state=mixscape_random_state)

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
    ms.perturbation_signature(
        adata,
        pert_key=pert_col,
        control=control_label,
        ref_selection_mode=signature_ref_selection_mode,
        split_by=signature_split_by,
        n_neighbors=signature_n_neighbors,
        use_rep="X_pca",
        n_dims=min(signature_n_dims, max_pcs),
        batch_size=batch_size,
    )
    mixscape_kwargs = dict(
        adata=adata,
        control=control_label,
        layer="X_pert",
        min_de_genes=mixscape_min_de_genes,
        logfc_threshold=mixscape_logfc_threshold,
        test_method=mixscape_test_method,
        iter_num=mixscape_iter_num,
        scale=mixscape_scale,
        split_by=mixscape_split_by,
        pval_cutoff=mixscape_pval_cutoff,
        perturbation_type=mixscape_perturbation_type,
        random_state=mixscape_random_state,
    )
    if "pert_key" in inspect.signature(ms.mixscape).parameters:
        mixscape_kwargs["pert_key"] = pert_col
    else:
        mixscape_kwargs["labels"] = pert_col
    ms.mixscape(**mixscape_kwargs)
    return {"pca_n_comps": max_pcs, "signature_n_dims": min(signature_n_dims, max_pcs), "batch_size": batch_size}


def summarize_chunk(
    adata: ad.AnnData,
    pert_col: str,
    control_label: str,
    pca_dims: int,
    perturbed_class: str = "KO",
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
    perturbed_mask = detect_perturbed_mask(obs, perturbed_class)

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
        pert_col=args.pert_col,
        csc_max_genes=args.csc_max_genes,
        csc_max_total_nnz=args.csc_max_total_nnz,
    )
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")

    batch_size = args.batch_size if args.batch_size and args.batch_size > 0 else None

    effective_parameters = run_mixscape(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        use_hvg_for_pca=args.use_hvg_for_pca,
        pca_n_comps=args.pca_n_comps,
        signature_n_neighbors=args.signature_n_neighbors,
        signature_n_dims=args.signature_n_dims,
        signature_ref_selection_mode=args.signature_ref_selection_mode,
        signature_split_by=args.signature_split_by or None,
        batch_size=batch_size,
        auto_batch_max_elements=args.auto_batch_max_elements,
        auto_batch_size=args.auto_batch_size,
        normalize_target_sum=args.normalize_target_sum,
        mixscape_logfc_threshold=args.mixscape_logfc_threshold,
        mixscape_pval_cutoff=args.mixscape_pval_cutoff,
        mixscape_min_de_genes=args.mixscape_min_de_genes,
        mixscape_test_method=args.mixscape_test_method,
        mixscape_iter_num=args.mixscape_iter_num,
        mixscape_scale=args.mixscape_scale,
        mixscape_split_by=args.mixscape_split_by or None,
        mixscape_perturbation_type=args.mixscape_perturbation_type,
        mixscape_random_state=args.mixscape_random_state,
    )
    stats_df, effects_df, labels_df = summarize_chunk(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        pca_dims=args.pca_dims,
        perturbed_class=args.mixscape_perturbation_type,
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
        **runtime_metadata(adata),
        "effective_parameters": effective_parameters,
        "effect_coordinate_system": "chunk-specific PCA; do not compare components across chunks",
        "h5ad": str(args.h5ad),
        "chunk_cells": str(args.chunk_cells),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "pert_col": args.pert_col,
        "control_label": args.control_label,
        "pca_dims": args.pca_dims,
        "pca_n_comps": args.pca_n_comps,
        "signature_n_neighbors": args.signature_n_neighbors,
        "signature_n_dims": args.signature_n_dims,
        "signature_ref_selection_mode": args.signature_ref_selection_mode,
        "signature_split_by": args.signature_split_by or None,
        "use_hvg_for_pca": bool(args.use_hvg_for_pca),
        "batch_size": batch_size,
        "auto_batch_max_elements": int(args.auto_batch_max_elements),
        "auto_batch_size": int(args.auto_batch_size),
        "csc_max_genes": int(args.csc_max_genes),
        "csc_max_total_nnz": int(args.csc_max_total_nnz),
        "normalize_target_sum": args.normalize_target_sum,
        "mixscape_logfc_threshold": args.mixscape_logfc_threshold,
        "mixscape_pval_cutoff": args.mixscape_pval_cutoff,
        "mixscape_min_de_genes": args.mixscape_min_de_genes,
        "mixscape_test_method": args.mixscape_test_method,
        "mixscape_iter_num": args.mixscape_iter_num,
        "mixscape_scale": bool(args.mixscape_scale),
        "mixscape_split_by": args.mixscape_split_by or None,
        "mixscape_perturbation_type": args.mixscape_perturbation_type,
        "mixscape_random_state": args.mixscape_random_state,
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
    ap.add_argument("--pca-n-comps", type=int, default=60)
    ap.add_argument("--signature-n-neighbors", type=int, default=30)
    ap.add_argument("--signature-n-dims", type=int, default=60)
    ap.add_argument("--signature-ref-selection-mode", choices=("nn", "split_by"), default="nn")
    ap.add_argument("--signature-split-by", default="")
    ap.add_argument("--mixscape-min-de-genes", type=int, default=1)
    ap.add_argument("--mixscape-logfc-threshold", type=float, default=0.10)
    ap.add_argument("--mixscape-pval-cutoff", type=float, default=0.05)
    ap.add_argument("--mixscape-test-method", choices=("logreg", "t-test", "wilcoxon", "t-test_overestim_var"), default="wilcoxon")
    ap.add_argument("--mixscape-iter-num", type=int, default=10)
    ap.add_argument("--mixscape-scale", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--mixscape-split-by", default="")
    ap.add_argument("--mixscape-perturbation-type", default="KO")
    ap.add_argument("--mixscape-random-state", type=int, default=0)
    return ap.parse_args()

def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        h5ad=Path(str(snk.input.h5ad)),
        chunk_cells=Path(str(snk.input.chunk_cells)),
        output_dir=Path(str(snk.params.outdir)),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        pca_dims=int(snk.params.pca_dims),
        pca_n_comps=int(snk.params.pca_n_comps),
        signature_n_neighbors=int(snk.params.signature_n_neighbors),
        signature_n_dims=int(snk.params.signature_n_dims),
        signature_ref_selection_mode=str(snk.params.signature_ref_selection_mode),
        signature_split_by=str(snk.params.signature_split_by),
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
        mixscape_min_de_genes=int(snk.params.min_de_genes),
        mixscape_test_method=str(snk.params.test_method),
        mixscape_iter_num=int(snk.params.iter_num),
        mixscape_scale=bool(snk.params.scale),
        mixscape_split_by=str(snk.params.split_by),
        mixscape_perturbation_type=str(snk.params.perturbation_type),
        mixscape_random_state=int(snk.params.random_state),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
