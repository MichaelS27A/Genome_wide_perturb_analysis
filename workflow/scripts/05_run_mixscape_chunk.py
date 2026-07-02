#!/usr/bin/env python3
"""Run chunk-level Mixscape and produce perturbation effect summaries."""

from __future__ import annotations

import argparse
import inspect
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
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


def load_subset(h5ad_path: Path, chunk_cells_tsv: Path) -> ad.AnnData:
    chunk_df = pd.read_csv(chunk_cells_tsv, sep="\t", compression="infer", dtype="string")
    barcodes = chunk_df["cell_barcode"].dropna().astype(str).unique().tolist()

    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs_names = pd.Index(adata.obs_names.astype(str))
    keep = obs_names.isin(barcodes)
    if int(keep.sum()) == 0:
        adata.file.close()
        raise RuntimeError("No overlap between chunk barcodes and adata.obs_names")

    try:
        sub = adata[keep].to_memory()
    except Exception as e:
        idx = np.where(keep)[0]
        batch = 512
        parts: list[ad.AnnData] = []
        for i in range(0, len(idx), batch):
            j = min(i + batch, len(idx))
            parts.append(adata[idx[i:j]].to_memory())
        if not parts:
            raise RuntimeError("Subset fallback failed: no cells loaded") from e
        sub = ad.concat(parts, axis=0, join="outer", merge="same")
    finally:
        adata.file.close()

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


def main() -> None:
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
    ap.add_argument("--write-subset", action="store_true")
    ap.add_argument(
        "--use-hvg-for-pca",
        action="store_true",
        help="If set, PCA is computed on highly variable genes only. Default uses full transcriptome.",
    )
    ap.add_argument("--normalize-target-sum", type=float, default=1e4)
    ap.add_argument("--mixscape-logfc-threshold", type=float, default=0.10)
    ap.add_argument("--mixscape-pval-cutoff", type=float, default=0.05)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_subset(args.h5ad, args.chunk_cells)
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
        "normalize_target_sum": args.normalize_target_sum,
        "mixscape_logfc_threshold": args.mixscape_logfc_threshold,
        "mixscape_pval_cutoff": args.mixscape_pval_cutoff,
        "n_selected_perturbed_cells": int(sel.shape[0]),
    }
    (args.output_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
    (args.output_dir / "done.txt").write_text("ok\n")

    print(f"Wrote {args.output_dir / 'done.txt'}")


if __name__ == "__main__":
    main()
