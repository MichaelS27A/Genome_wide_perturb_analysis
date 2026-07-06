#!/usr/bin/env python3
"""X-Atlas sgRNA-abundance stratification diagnostic against Mixscape calls.

This script reproduces the key preprint-style idea (low vs high sgRNA abundance
within each perturbation) using available local fields:
  - per-cell guide abundance from guide_calls CSV (`num_umis`)
  - per-cell KO labels from finished Mixscape chunk outputs
"""


import argparse
import json
import glob
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact, spearmanr


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labels-glob", type=str, required=True, help="Glob for cell_mixscape_labels.tsv.gz files")
    ap.add_argument("--guide-calls-csv", type=Path, required=True, help="Guide-calls CSV(.gz) path")
    ap.add_argument("--outdir", type=Path, required=True, help="Output directory")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--min-cells-per-perturbation", type=int, default=100)
    ap.add_argument("--min-cells-per-stratum", type=int, default=20)
    ap.add_argument("--guide-chunksize", type=int, default=400_000)
    return ap.parse_args()


def detect_ko(obs: pd.DataFrame) -> pd.Series:
    if "is_predicted_perturbed" in obs.columns:
        ko = pd.to_numeric(obs["is_predicted_perturbed"], errors="coerce").fillna(0).astype(int) > 0
    else:
        ko = pd.Series(False, index=obs.index)

    if "mixscape_class_global" in obs.columns:
        ko2 = obs["mixscape_class_global"].astype(str).str.lower().str.contains("ko|perturbed", regex=True, na=False)
        ko = ko | ko2
    elif "mixscape_class" in obs.columns:
        ko2 = obs["mixscape_class"].astype(str).str.lower().str.contains("ko|perturbed", regex=True, na=False)
        ko = ko | ko2

    return ko.astype(int)


def load_mixscape_labels(labels_glob: str) -> pd.DataFrame:
    files = [Path(p) for p in sorted(glob.glob(labels_glob))]
    if not files:
        raise FileNotFoundError(f"No files matched: {labels_glob}")

    dfs = []
    for p in files:
        done_file = p.parent / "done.txt"
        if not done_file.exists():
            continue
        df = pd.read_csv(
            p,
            sep="\t",
            compression="infer",
            usecols=lambda c: c in {
                "cell_barcode",
                "gene_target",
                "is_predicted_perturbed",
                "mixscape_class_global",
                "mixscape_class",
            },
            dtype="string",
        )
        df["source_chunk"] = p.parent.name
        dfs.append(df)

    if not dfs:
        raise RuntimeError("No finished chunk label files found (requires done.txt in chunk directory).")

    labels = pd.concat(dfs, axis=0, ignore_index=True)
    labels["cell_barcode"] = labels["cell_barcode"].astype(str)
    labels["gene_target"] = labels["gene_target"].astype(str)
    labels["is_ko"] = detect_ko(labels)

    # Drop duplicated barcodes (controls are intentionally repeated across chunks).
    labels = labels.sort_values("source_chunk").drop_duplicates(subset=["cell_barcode"], keep="first")
    return labels[["cell_barcode", "gene_target", "is_ko", "source_chunk"]].copy()


def total_guide_umi(series: pd.Series) -> pd.Series:
    # `num_umis` stores dual-guide counts like "1794|1481".
    # Sum all integer tokens found in each cell.
    extracted = series.astype(str).str.extractall(r"(\d+)")
    if extracted.empty:
        return pd.Series(0, index=series.index, dtype=np.int64)
    vals = extracted[0].astype(np.int64).groupby(level=0).sum()
    return vals.reindex(series.index).fillna(0).astype(np.int64)


def load_filtered_guide_calls(guide_csv: Path, keep_barcodes: set[str], chunksize: int) -> pd.DataFrame:
    usecols = ["cell_barcode", "gene_target", "num_umis", "pass_guide_filter"]
    parts = []

    for chunk in pd.read_csv(
        guide_csv,
        compression="infer",
        usecols=usecols,
        dtype={"cell_barcode": "string", "gene_target": "string", "num_umis": "string"},
        chunksize=chunksize,
    ):
        chunk["cell_barcode"] = chunk["cell_barcode"].astype(str)
        chunk = chunk[chunk["cell_barcode"].isin(keep_barcodes)]
        if chunk.empty:
            continue
        if "pass_guide_filter" in chunk.columns:
            chunk = chunk[chunk["pass_guide_filter"].astype(str).str.lower().eq("true")]
        if chunk.empty:
            continue

        chunk["total_guide_umis"] = total_guide_umi(chunk["num_umis"])
        parts.append(chunk[["cell_barcode", "gene_target", "total_guide_umis"]].copy())

    if not parts:
        return pd.DataFrame(columns=["cell_barcode", "gene_target", "total_guide_umis"])

    out = pd.concat(parts, axis=0, ignore_index=True)
    out = out.drop_duplicates(subset=["cell_barcode"], keep="first")
    out["total_guide_umis"] = pd.to_numeric(out["total_guide_umis"], errors="coerce").fillna(0).astype(np.int64)
    return out


def per_perturbation_stratification(
    df: pd.DataFrame,
    min_cells: int,
    min_per_group: int,
    control_label: str,
) -> pd.DataFrame:
    rows = []

    nonctrl = df[df["gene_target"] != control_label].copy()
    for pert, g in nonctrl.groupby("gene_target", sort=True):
        n = len(g)
        if n < min_cells:
            continue

        q25 = float(np.quantile(g["total_guide_umis"], 0.25))
        q75 = float(np.quantile(g["total_guide_umis"], 0.75))
        low = g[g["total_guide_umis"] <= q25]
        high = g[g["total_guide_umis"] >= q75]

        n_low = int(len(low))
        n_high = int(len(high))
        if n_low < min_per_group or n_high < min_per_group:
            continue

        ko_low = int(low["is_ko"].sum())
        ko_high = int(high["is_ko"].sum())
        rate_low = ko_low / n_low if n_low else np.nan
        rate_high = ko_high / n_high if n_high else np.nan

        table = np.array([[ko_high, n_high - ko_high], [ko_low, n_low - ko_low]], dtype=int)
        try:
            odds_ratio, p_greater = fisher_exact(table, alternative="greater")
        except Exception:
            odds_ratio, p_greater = np.nan, np.nan

        rows.append(
            {
                "gene_target": pert,
                "n_cells": n,
                "q25_guide_umi": q25,
                "q75_guide_umi": q75,
                "n_low": n_low,
                "n_high": n_high,
                "ko_low": ko_low,
                "ko_high": ko_high,
                "ko_rate_low": rate_low,
                "ko_rate_high": rate_high,
                "ko_rate_delta_high_minus_low": rate_high - rate_low,
                "ko_rate_ratio_high_over_low": (rate_high / rate_low) if rate_low > 0 else np.nan,
                "median_guide_umi_low": float(np.median(low["total_guide_umis"])) if n_low else np.nan,
                "median_guide_umi_high": float(np.median(high["total_guide_umis"])) if n_high else np.nan,
                "fisher_odds_ratio_high_gt_low": odds_ratio,
                "fisher_pvalue_high_gt_low": p_greater,
            }
        )

    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("ko_rate_delta_high_minus_low", ascending=False).reset_index(drop=True)
    return out


def write_summary(
    joined: pd.DataFrame,
    pert_df: pd.DataFrame,
    outdir: Path,
    control_label: str,
) -> None:
    nonctrl = joined[joined["gene_target"] != control_label].copy()
    rho, rho_p = spearmanr(np.log1p(nonctrl["total_guide_umis"]), nonctrl["is_ko"]) if not nonctrl.empty else (np.nan, np.nan)

    dec = nonctrl.copy()
    dec["umi_decile"] = pd.qcut(dec["total_guide_umis"], 10, labels=False, duplicates="drop")
    decile = (
        dec.groupby("umi_decile", as_index=False)
        .agg(
            n_cells=("cell_barcode", "size"),
            n_ko=("is_ko", "sum"),
            median_guide_umi=("total_guide_umis", "median"),
            mean_guide_umi=("total_guide_umis", "mean"),
        )
        .sort_values("umi_decile")
    )
    decile["ko_rate"] = decile["n_ko"] / decile["n_cells"]
    decile.to_csv(outdir / "hct_mixscape_ko_by_umi_decile.tsv", sep="\t", index=False)

    summary = {
        "n_cells_joined_total": int(len(joined)),
        "n_cells_joined_noncontrol": int((joined["gene_target"] != control_label).sum()),
        "n_perturbations_tested": int(len(pert_df)),
        "weighted_ko_rate_low": float(pert_df["ko_low"].sum() / pert_df["n_low"].sum()) if len(pert_df) else np.nan,
        "weighted_ko_rate_high": float(pert_df["ko_high"].sum() / pert_df["n_high"].sum()) if len(pert_df) else np.nan,
        "fraction_perturbations_high_gt_low": float((pert_df["ko_rate_delta_high_minus_low"] > 0).mean()) if len(pert_df) else np.nan,
        "fraction_perturbations_delta_ge_0p1": float((pert_df["ko_rate_delta_high_minus_low"] >= 0.10).mean()) if len(pert_df) else np.nan,
        "median_delta_high_minus_low": float(pert_df["ko_rate_delta_high_minus_low"].median()) if len(pert_df) else np.nan,
        "spearman_log1p_umi_vs_is_ko": float(rho) if pd.notna(rho) else np.nan,
        "spearman_pvalue": float(rho_p) if pd.notna(rho_p) else np.nan,
    }
    (outdir / "hct_sgrna_stratified_mixscape_summary.json").write_text(json.dumps(summary, indent=2))


def make_plots(pert_df: pd.DataFrame, decile_df: pd.DataFrame, outdir: Path) -> None:
    sns.set_theme(style="whitegrid")

    # Paired KO rates (low vs high)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(pert_df["ko_rate_low"], pert_df["ko_rate_high"], s=10, alpha=0.4)
    ax.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("KO rate (low sgRNA quartile)")
    ax.set_ylabel("KO rate (high sgRNA quartile)")
    ax.set_title("HCT116: KO-rate shift with sgRNA abundance")
    fig.tight_layout()
    fig.savefig(outdir / "hct_ko_rate_low_vs_high_sgrna_quartiles.png", dpi=180)
    plt.close(fig)

    # Delta distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(pert_df["ko_rate_delta_high_minus_low"], bins=40, color="#2a6f97", alpha=0.9)
    ax.axvline(0.0, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("KO-rate delta (high - low sgRNA quartile)")
    ax.set_ylabel("Number of perturbations")
    ax.set_title("HCT116: Distribution of sgRNA-stratified KO-rate deltas")
    fig.tight_layout()
    fig.savefig(outdir / "hct_ko_rate_delta_high_minus_low_hist.png", dpi=180)
    plt.close(fig)

    # KO rate by global UMI decile
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(decile_df["umi_decile"], decile_df["ko_rate"], marker="o", linewidth=2)
    ax.set_xlabel("Global sgRNA UMI decile (non-control cells)")
    ax.set_ylabel("KO rate")
    ax.set_title("HCT116: Mixscape KO rate by sgRNA UMI decile")
    fig.tight_layout()
    fig.savefig(outdir / "hct_mixscape_ko_rate_by_umi_decile.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    labels = load_mixscape_labels(args.labels_glob)
    keep_barcodes = set(labels["cell_barcode"].tolist())
    guide = load_filtered_guide_calls(args.guide_calls_csv, keep_barcodes, args.guide_chunksize)

    joined = labels.merge(guide, on="cell_barcode", how="inner", suffixes=("_mixscape", "_guide"))
    if joined.empty:
        raise RuntimeError("No overlap between labels and guide calls after filtering.")

    # Prefer Mixscape perturbation identity; guide target retained for QA only.
    joined["gene_target"] = joined["gene_target_mixscape"].astype(str)
    joined["gene_target_guide"] = joined["gene_target_guide"].astype(str)
    joined["gene_target_match"] = joined["gene_target"] == joined["gene_target_guide"]

    qa = {
        "n_label_cells_after_dedup": int(len(labels)),
        "n_joined_cells": int(len(joined)),
        "fraction_gene_target_match_mixscape_vs_guide_calls": float(joined["gene_target_match"].mean()),
        "n_nonmatching_gene_target": int((~joined["gene_target_match"]).sum()),
    }
    (args.outdir / "hct_sgrna_join_qc.json").write_text(json.dumps(qa, indent=2))

    pert = per_perturbation_stratification(
        joined[["cell_barcode", "gene_target", "is_ko", "total_guide_umis"]].copy(),
        min_cells=args.min_cells_per_perturbation,
        min_per_group=args.min_cells_per_stratum,
        control_label=args.control_label,
    )
    pert.to_csv(args.outdir / "hct_sgrna_stratified_mixscape_per_perturbation.tsv", sep="\t", index=False)

    write_summary(joined, pert, args.outdir, args.control_label)
    decile_df = pd.read_csv(args.outdir / "hct_mixscape_ko_by_umi_decile.tsv", sep="\t")

    if not pert.empty and not decile_df.empty:
        make_plots(pert, decile_df, args.outdir)

    run_meta = {
        "labels_glob": args.labels_glob,
        "guide_calls_csv": str(args.guide_calls_csv),
        "outdir": str(args.outdir),
        "control_label": args.control_label,
        "min_cells_per_perturbation": args.min_cells_per_perturbation,
        "min_cells_per_stratum": args.min_cells_per_stratum,
        "guide_chunksize": args.guide_chunksize,
    }
    (args.outdir / "run_meta.json").write_text(json.dumps(run_meta, indent=2))
    print(f"Wrote outputs to {args.outdir}")


if __name__ == "__main__":
    main()
