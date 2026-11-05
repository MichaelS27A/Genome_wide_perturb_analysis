#!/usr/bin/env python3
"""Cluster complete saved DE vectors from historical Mixscape-selected cells."""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/anka_signature_matplotlib")
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import h5py
from anndata.io import read_elem
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location(
    "signature_plotting", Path(__file__).with_name("51_cluster_custom_mixscape_signatures.py")
)
plotting = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plotting)

DATASETS = ("STEMCELL", "HEK293T", "HCT116", "K562_GWPS_RAW_FULL")


def dataset_analysis(workspace, dataset, out):
    out.mkdir(parents=True, exist_ok=True)
    base = workspace / "X_Atlas/results/mixscape_pipeline" / dataset
    old = base / "custom_target_heatmaps/slc4_ampk_metabolism/mixscape"
    old_meta = json.loads((old / "run_meta.json").read_text())
    target_path = workspace / "X_Atlas/code/config/targets/slc4_ampk_metabolism_targets.txt"
    targets = target_path.read_text().splitlines()
    goi = set((target_path.parent / "slc4_ampk_metabolism_goi_targets.txt").read_text().splitlines())
    candidates = set(old_meta["selected_hypusine_perturbations"])
    old_randoms = set(old_meta["selected_random_controls"])
    files, available_targets, comparators = [], set(), set()
    for p in sorted((base / "postprocess/chunk_runs").glob("*/perturbation_differential_genes_full.tsv.gz")):
        chunk_targets = set((base / "chunks" / (p.parent.name + "_perturbations.txt")).read_text().splitlines())
        if chunk_targets & candidates:
            files.append(p)
            available_targets.update(chunk_targets & candidates)
            comparators.update(chunk_targets & old_randoms)
    wanted = available_targets | comparators
    counts = pd.read_csv(old / "assignment_counts.tsv", sep="\t").set_index("perturbation").n_assigned_cells
    parts = []
    source_rows = []
    columns = ["perturbation", "gene", "score", "logfoldchange", "pval_adj", "pct_nz_group"]
    for p in files:
        print(f"[saved DE] {dataset}: {p.parent.name}", flush=True)
        for frame in pd.read_csv(p, sep="\t", usecols=columns, chunksize=200000):
            keep = frame[frame.perturbation.isin(wanted)]
            if len(keep):
                parts.append(keep)
        source_rows.append({"path": str(p), "size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns})
    coverage = []
    for target in targets:
        n = int(counts.get(target, 0))
        status = (
            "no_mixscape_selected_cells" if n == 0 else
            "below_30_selected_cells" if n < 30 else
            "no_saved_full_DE_table" if target not in available_targets else "available"
        )
        coverage.append({"perturbation": target, "is_goi": target in goi, "n_selected_cells": n, "status": status})
    qc = pd.DataFrame(coverage)
    qc.to_csv(out / "target_coverage_qc.tsv", sep="\t", index=False)
    if not parts:
        summary = {"dataset": dataset, "status": "no_complete_saved_DE_profiles", "n_requested_profiles": 0}
        (out / "run_meta.json").write_text(json.dumps(summary, indent=2))
        return summary
    data = pd.concat(parts, ignore_index=True)
    if data.duplicated(["perturbation", "gene"]).any():
        raise ValueError("Duplicate target/readout pairs in complete DE tables")
    # A full table includes negative and nonsignificant effects. Absence is never zero-filled.
    lfc = data.pivot(index="perturbation", columns="gene", values="logfoldchange")
    score = data.pivot(index="perturbation", columns="gene", values="score").reindex_like(lfc)
    detection = data.pivot(index="perturbation", columns="gene", values="pct_nz_group").reindex_like(lfc)
    fdr = data.pivot(index="perturbation", columns="gene", values="pval_adj").reindex_like(lfc)
    complete = np.isfinite(lfc).all(axis=0) & np.isfinite(score).all(axis=0)
    # Match symbols to feature IDs, particularly the Ensembl-indexed K562 atlas.
    symbol = pd.Series(lfc.columns, index=lfc.columns)
    if dataset == "K562_GWPS_RAW_FULL":
        with h5py.File(Path(old_meta["h5ad"]).resolve(), "r") as h5:
            var = read_elem(h5["var"])
        symbol = var.gene_name.reindex(lfc.columns).fillna(symbol).astype(str)
    excluded = symbol.isin(set(targets) | comparators)
    # Restrict clustering to a coherent response panel, then retain every signed
    # value for those readouts, including values nonsignificant in other targets.
    regs = lfc.index.intersection(sorted(available_targets))
    response = (
        (fdr.loc[regs] < .05)
        & (lfc.loc[regs].abs() >= .25)
        & (detection.loc[regs] >= .10)
    ).any(axis=0)
    eligible = complete & response & ~excluded
    selected_features = score.loc[regs, eligible].abs().max(axis=0).nlargest(2000).index
    if len(selected_features) < 3:
        raise ValueError(f"{dataset}: fewer than three eligible response features")
    feature_qc = pd.DataFrame({
        "gene_id": lfc.columns, "gene_symbol": symbol.to_numpy(),
        "complete_values": complete.to_numpy(), "direct_target_excluded": excluded.to_numpy(),
        "significant_response_in_requested_target": response.to_numpy(),
        "used_for_clustering": lfc.columns.isin(selected_features),
    })
    feature_qc.to_csv(out / "feature_manifest.tsv", sep="\t", index=False)
    # Clipping is a recorded sensitivity choice to limit unstable huge ratios
    # when reference expression approaches zero; unmodified vectors are exported.
    lfc.to_csv(out / "full_unclipped_log2_foldchange_profiles.tsv.gz", sep="\t")
    score.to_csv(out / "full_signed_wilcoxon_profiles.tsv.gz", sep="\t")
    plot_labels = {
        g: str(symbol[g]) if not symbol.duplicated(keep=False)[g] else f"{symbol[g]} [{g}]"
        for g in selected_features
    }
    profile_qc = pd.DataFrame({
        "perturbation": lfc.index, "n_profile_cells": [int(counts.get(t, 0)) for t in lfc.index],
        "category": ["requested" if t in available_targets else "random_comparator" for t in lfc.index],
    })
    profile_qc["is_goi"] = profile_qc.perturbation.isin(goi)
    results = {}
    for name, frame in [
        ("logfoldchange_clustering", lfc.loc[:, selected_features].clip(-3, 3)),
        ("signed_score_sensitivity", score.loc[:, selected_features]),
    ]:
        folder = out / name
        folder.mkdir(exist_ok=True)
        frame = frame.rename(columns=plot_labels)
        results[name] = plotting.save_clustering(
            frame, profile_qc, folder,
            f"{dataset}: saved Mixscape-selected expression changes",
            effect_label="Log2 fold change versus NTC (clipped ±3)" if name.startswith("log") else "Signed Wilcoxon score versus NTC",
        )
    target_folder = out / "requested_targets_only"
    target_folder.mkdir(exist_ok=True)
    results["requested_targets_only"] = plotting.save_clustering(
        lfc.loc[regs, selected_features].clip(-3, 3).rename(columns=plot_labels),
        profile_qc, target_folder, f"{dataset}: requested targets, expression-change profiles",
        effect_label="Log2 fold change versus NTC (clipped ±3)",
    )
    broad_features = lfc.columns[complete & (detection >= .10).any(axis=0) & ~excluded]
    broad_folder = out / "broad_detected_feature_sensitivity"
    broad_folder.mkdir(exist_ok=True)
    results["broad_detected_feature_sensitivity"] = plotting.save_clustering(
        lfc.loc[:, broad_features].clip(-3, 3), profile_qc, broad_folder,
        f"{dataset}: all detected readouts, expression-change profiles",
        effect_label="Log2 fold change versus NTC (clipped ±3)",
    )
    # Compare the old expression and new effect correlations on the same target
    # pairs, rather than comparing different populations or feature universes.
    old_corr = pd.read_csv(old / "correlation_matrix.tsv.gz", sep="\t", index_col=0)
    absolute = pd.read_csv(old / "pseudobulk_log1p_matrix.tsv.gz", sep="\t", index_col=0)
    if not set(selected_features).issubset(absolute.columns):
        raise ValueError("Saved DE feature IDs do not match historical expression profiles")
    matched_absolute = absolute.loc[lfc.index, selected_features]
    matched_corr = pd.DataFrame(
        np.corrcoef(matched_absolute), index=lfc.index, columns=lfc.index,
    )
    matched_corr.to_csv(out / "matched_feature_absolute_expression_correlations.tsv", sep="\t")
    new_corr = pd.read_csv(out / "logfoldchange_clustering/profile_correlation_matrix.tsv", sep="\t", index_col=0)
    comparisons = []
    same_regs = sorted(set(regs) & set(new_corr.index) & set(old_corr.index))
    for i, first in enumerate(same_regs):
        for second in same_regs[i+1:]:
            comparisons.append({
                "target_a": first, "target_b": second,
                "original_absolute_expression_correlation": old_corr.loc[first, second],
                "absolute_expression_correlation_on_same_readouts": matched_corr.loc[first, second],
                "saved_effect_profile_correlation": new_corr.loc[first, second],
            })
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(out / "matched_requested_pair_comparison.tsv", sep="\t", index=False)
    qc["status"] = np.where(qc.perturbation.isin(regs), "clustered", qc.status)
    qc.to_csv(out / "target_coverage_qc.tsv", sep="\t", index=False)
    summary = {
        "dataset": dataset, "status": "complete", "n_requested_profiles": len(regs),
        "n_comparator_profiles": len(lfc) - len(regs), "n_clustering_features": len(selected_features),
        "requested_targets": list(regs), "comparators": sorted(comparators),
        "method": "saved full gene-level DE from historical Mixscape-selected cells versus NTCs",
        "is_native_X_pert": False, "recomputed_signatures_or_classifications": False,
        "linkage": "average linkage on Pearson correlation distance",
        "feature_selection": "up to 2000 highest-|Wilcoxon score| readouts with FDR<0.05, |log2FC|>=0.25 and target detection>=10% in at least one requested target",
        "direct_targets_excluded": sorted(set(targets) | comparators),
        "profile_values": "all signed log2FC values on shared selected readouts; clipped to [-3,3]",
        "sensitivity_values": "all signed Wilcoxon scores on the same readouts",
        "missing_values": "features with any nonfinite/missing effect removed; never zero-imputed",
        "comparators_note": "previously sampled random perturbations with full profiles in the same processing chunks as requested targets",
        "sources": source_rows, "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "results": results,
        "median_original_requested_pair_correlation": float(comparison.original_absolute_expression_correlation.median()) if len(comparison) else None,
        "median_matched_readout_absolute_expression_correlation": float(comparison.absolute_expression_correlation_on_same_readouts.median()) if len(comparison) else None,
        "median_saved_effect_requested_pair_correlation": float(comparison.saved_effect_profile_correlation.median()) if len(comparison) else None,
        "limitations": "Descriptive, selected-cell analysis; feature selection is post hoc. Clusters are not evidence of pathway specificity or biological replication. Full profiles are unavailable for some targets.",
    }
    (out / "run_meta.json").write_text(json.dumps(summary, indent=2))
    print("[complete]", dataset, len(regs), "requested profiles", len(selected_features), "features", flush=True)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    summaries = [dataset_analysis(args.workspace, ds, args.outdir / ds) for ds in DATASETS]
    (args.outdir / "summary.json").write_text(json.dumps(summaries, indent=2))
    pd.DataFrame([{
        key: row.get(key) for key in (
            "dataset", "status", "n_requested_profiles", "n_comparator_profiles", "n_clustering_features",
            "median_original_requested_pair_correlation", "median_saved_effect_requested_pair_correlation",
            "median_matched_readout_absolute_expression_correlation",
        )
    } for row in summaries]).to_csv(args.outdir / "dataset_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
