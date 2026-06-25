#!/usr/bin/env python3
"""Run chunk-level Mixscape and produce perturbation effect summaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats


def detect_perturbed_mask(obs: pd.DataFrame) -> pd.Series:
    """Heuristic for perturbed-class cells across different Mixscape label variants."""
    candidates = [c for c in obs.columns if "mixscape" in c.lower() and "class" in c.lower()]
    if not candidates:
        return pd.Series(False, index=obs.index)

    col = candidates[0]
    vals = obs[col].astype(str).str.lower()
    # Common conventions: KO vs NP, or 'perturbed' strings.
    mask = vals.str.contains("ko|perturbed", regex=True, na=False)
    return mask


def load_subset(h5ad_path: Path, chunk_cells_tsv: Path) -> ad.AnnData:
    chunk_df = pd.read_csv(chunk_cells_tsv, sep="\t", compression="infer", dtype="string")
    barcodes = chunk_df["cell_barcode"].dropna().astype(str).unique().tolist()

    adata = sc.read_h5ad(h5ad_path, backed="r")
    obs_names = pd.Index(adata.obs_names.astype(str))
    keep = obs_names.isin(barcodes)
    if int(keep.sum()) == 0:
        adata.file.close()
        raise RuntimeError("No overlap between chunk barcodes and adata.obs_names")

    sub = adata[keep].to_memory()
    adata.file.close()

    chunk_meta = chunk_df.drop_duplicates(subset=["cell_barcode"]).set_index("cell_barcode")
    common = sub.obs_names.intersection(chunk_meta.index)
    sub.obs = sub.obs.copy()
    sub.obs.loc[common, "chunk_gene_target"] = chunk_meta.loc[common, "gene_target"].astype(str).values
    sub.obs.loc[common, "chunk_is_control"] = chunk_meta.loc[common, "is_control"].astype(int).values
    return sub


def run_mixscape(
    adata: ad.AnnData, pert_col: str, control_label: str, use_hvg_for_pca: bool = False
) -> None:
    import pertpy as pt

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # Keep HVG annotation for diagnostics, but default to full-transcriptome PCA.
    sc.pp.highly_variable_genes(adata, subset=False)
    sc.pp.pca(adata, use_highly_variable=use_hvg_for_pca)

    ms = pt.tl.Mixscape()
    split_by = "sample" if "sample" in adata.obs.columns else None
    ms.perturbation_signature(adata, pert_key=pert_col, control=control_label, split_by=split_by)
    ms.mixscape(adata=adata, control=control_label, labels=pert_col, layer="X_pert")


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

    # Cell-level output with perturbed flags.
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

        # Per-PC one-sample difference vs control via two-sample t-test.
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
    ap.add_argument("--write-subset", action="store_true")
    ap.add_argument(
        "--use-hvg-for-pca",
        action="store_true",
        help="If set, PCA is computed on highly variable genes only. Default uses full transcriptome.",
    )
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    adata = load_subset(args.h5ad, args.chunk_cells)
    if args.pert_col not in adata.obs.columns:
        raise RuntimeError(f"perturbation column '{args.pert_col}' not found in adata.obs")

    run_mixscape(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        use_hvg_for_pca=args.use_hvg_for_pca,
    )
    stats_df, effects_df, labels_df = summarize_chunk(
        adata,
        pert_col=args.pert_col,
        control_label=args.control_label,
        pca_dims=args.pca_dims,
    )

    # Keep only predicted-perturbed cells from non-control perturbations.
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
        "n_selected_perturbed_cells": int(sel.shape[0]),
    }
    (args.output_dir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
    (args.output_dir / "done.txt").write_text("ok\n")

    print(f"Wrote {args.output_dir / 'done.txt'}")


if __name__ == "__main__":
    main()
