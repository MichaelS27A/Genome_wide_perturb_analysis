#!/usr/bin/env python3
"""Recompute shared-reference Mixscape signatures and cluster target profiles.

Historical KO labels are retained. This does not reproduce chunk-specific PCA
or rerun classification. Each dataset has its own gene coordinates and clusters.
"""
import os
for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(key, "4")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/anka_signature_matplotlib")

import argparse
import hashlib
import json
import time
from pathlib import Path

import anndata as ad
from anndata.io import read_elem
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pertpy as pt
import scanpy as sc
from scipy import sparse
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform

DATASETS = ("STEMCELL", "HEK293T", "HCT116", "K562_GWPS_RAW_FULL")


def read_rows(h5, rows, genes):
    """Read selected raw-count rows without loading the full atlas matrix."""
    x = h5["X"]
    if isinstance(x, h5py.Dataset):
        return sparse.csr_matrix(x[rows, :][:, genes], dtype=np.float32)
    encoding = x.attrs.get("encoding-type", "")
    if encoding == "csc_matrix":
        mapping = np.full(int(x.attrs["shape"][0]), -1, dtype=np.int32)
        mapping[rows] = np.arange(len(rows))
        ptr = x["indptr"][:]
        data, indices, outptr = [], [], [0]
        for n, g in enumerate(genes):
            lo, hi = int(ptr[g]), int(ptr[g + 1])
            idx = mapping[x["indices"][lo:hi]]
            keep = idx >= 0
            values = x["data"][lo:hi][keep].astype(np.float32)
            indices.append(idx[keep])
            data.append(values)
            outptr.append(outptr[-1] + len(values))
            if n % 500 == 0:
                print(f"[read CSC] {n}/{len(genes)} genes", flush=True)
        return sparse.csc_matrix(
            (np.concatenate(data), np.concatenate(indices), np.asarray(outptr)),
            shape=(len(rows), len(genes)),
        ).tocsr()
    if encoding != "csr_matrix":
        raise ValueError(f"Unsupported matrix encoding {encoding!r}")
    ptr = x["indptr"][:]
    mapping = np.full(int(x.attrs["shape"][1]), -1, dtype=np.int32)
    mapping[genes] = np.arange(len(genes))
    data, indices, outptr = [], [], [0]
    for n, row in enumerate(rows):
        lo, hi = int(ptr[row]), int(ptr[row + 1])
        idx = mapping[x["indices"][lo:hi]]
        keep = idx >= 0
        values = x["data"][lo:hi][keep].astype(np.float32)
        indices.append(idx[keep])
        data.append(values)
        outptr.append(outptr[-1] + len(values))
        if n % 5000 == 0:
            print(f"[read CSR] {n}/{len(rows)} cells", flush=True)
    return sparse.csr_matrix(
        (np.concatenate(data), np.concatenate(indices), np.asarray(outptr)),
        shape=(len(rows), len(genes)),
    )


def save_clustering(profiles, qc, out, title, effect_label="Mean log1p expression residual versus nearest controls"):
    """Average linkage on correlation distance; omit undefined profiles."""
    valid = np.isfinite(profiles).all(axis=1) & (profiles.std(axis=1) > 1e-8)
    profiles = profiles.loc[valid]
    profiles.to_csv(out / "effect_profile_matrix.tsv.gz", sep="\t")
    if len(profiles) < 3:
        (out / "status.json").write_text(json.dumps({"status": "fewer than three valid profiles"}))
        return None
    distances = pdist(profiles.to_numpy(), metric="correlation")
    tree = linkage(distances, method="average", optimal_ordering=True)
    order = leaves_list(tree)
    corr = 1 - squareform(distances)
    np.fill_diagonal(corr, 1)
    correlation = pd.DataFrame(corr, index=profiles.index, columns=profiles.index)
    correlation.to_csv(out / "profile_correlation_matrix.tsv", sep="\t")
    pd.DataFrame(tree, columns=["left", "right", "distance", "n_members"]).to_csv(
        out / "linkage.tsv", sep="\t", index=False
    )
    assignments = qc.set_index("perturbation").loc[profiles.index].copy()
    assignments["leaf_order"] = np.argsort(order)
    for k in (2, 3, 4, 5, 6):
        if k <= len(profiles):
            assignments[f"cluster_k{k}"] = fcluster(tree, k, criterion="maxclust")
    assignments.sort_values("leaf_order").to_csv(out / "cluster_assignments.tsv", sep="\t")
    labels = [
        f"{p} ({int(assignments.loc[p, 'n_profile_cells'])})"
        for p in profiles.index[order]
    ]
    fig = plt.figure(figsize=(max(8, len(labels) * .17), max(7, len(labels) * .17)))
    grid = fig.add_gridspec(2, 2, width_ratios=[20, 1], height_ratios=[1, 5], hspace=.02)
    da = fig.add_subplot(grid[0, 0])
    dendrogram(tree, ax=da, no_labels=True, color_threshold=0, above_threshold_color="#475569")
    da.axis("off")
    ax = fig.add_subplot(grid[1, 0])
    im = ax.imshow(corr[np.ix_(order, order)], vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(np.arange(len(labels)), labels, rotation=90, fontsize=6)
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=6)
    colors = [
        "#1d4ed8" if bool(assignments.loc[p].get("is_goi", False)) else
        "#92400e" if assignments.loc[p, "category"] == "requested" else "#64748b"
        for p in profiles.index[order]
    ]
    for tick, color in zip(ax.get_xticklabels(), colors):
        tick.set_color(color)
    for tick, color in zip(ax.get_yticklabels(), colors):
        tick.set_color(color)
    fig.colorbar(im, cax=fig.add_subplot(grid[1, 1]), label="Pearson r of effect vectors")
    fig.suptitle(title, fontsize=12)
    fig.text(.23, .02, "Blue: AMPK/SLC4A4 targets; brown: other requested targets; grey: random comparators. Parentheses: selected-cell count.", fontsize=7)
    fig.subplots_adjust(bottom=.23, left=.23, top=.9)
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"clustered_effect_correlations.{suffix}", dpi=180)
    plt.close(fig)
    # Show signed effects directly, using the same perturbation order.
    variance = profiles.var(axis=0)
    display_genes = variance.nlargest(min(60, profiles.shape[1])).index
    display = profiles.loc[profiles.index[order], display_genes]
    gene_order = leaves_list(linkage(pdist(display.T, "euclidean"), "average"))
    display = display.iloc[:, gene_order]
    limit = max(.05, float(np.quantile(np.abs(display.to_numpy()), .98)))
    fig, ax = plt.subplots(figsize=(13, max(5, len(labels) * .18)))
    im = ax.imshow(display, aspect="auto", vmin=-limit, vmax=limit, cmap="RdBu_r")
    ax.set_yticks(np.arange(len(labels)), labels, fontsize=6)
    ax.set_xticks(np.arange(len(display.columns)), display.columns, rotation=90, fontsize=6)
    ax.set_title(title + "\nGenes with most variable effect profiles")
    fig.colorbar(im, ax=ax, label=effect_label)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"clustered_signed_effects.{suffix}", dpi=180)
    plt.close(fig)
    upper = np.triu_indices(len(profiles), 1)
    categories = assignments["category"].to_numpy()
    target_pair = (categories[upper[0]] == "requested") & (categories[upper[1]] == "requested")
    random_pair = (categories[upper[0]] == "random_comparator") & (categories[upper[1]] == "random_comparator")
    return {
        "n_profiles": len(profiles),
        "n_features": profiles.shape[1],
        "median_pairwise_effect_correlation": float(np.median(corr[upper])),
        "median_requested_pair_correlation": float(np.median(corr[upper][target_pair])) if target_pair.any() else None,
        "median_random_pair_correlation": float(np.median(corr[upper][random_pair])) if random_pair.any() else None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=DATASETS, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--max-csc-features", type=int, default=6000)
    parser.add_argument("--random-comparators", type=int, default=50)
    parser.add_argument("--min-selected-cells", type=int, default=30)
    args = parser.parse_args()
    started = time.time()
    args.outdir.mkdir(parents=True, exist_ok=True)
    base = args.workspace / "X_Atlas/results/mixscape_pipeline" / args.dataset
    target_file = args.workspace / "X_Atlas/code/config/targets/slc4_ampk_metabolism_targets.txt"
    targets = target_file.read_text().splitlines()
    goi = set((target_file.parent / "slc4_ampk_metabolism_goi_targets.txt").read_text().splitlines())
    old = base / "custom_target_heatmaps/slc4_ampk_metabolism/mixscape"
    old_meta = json.loads((old / "run_meta.json").read_text())
    rng = np.random.default_rng(0)
    randoms = sorted(rng.choice(old_meta["selected_random_controls"], args.random_comparators, replace=False))
    requested = set(targets) | set(randoms)
    chunk_rows = []
    for p in sorted((base / "chunks").glob("*_perturbations.txt")):
        if requested & set(p.read_text().splitlines()):
            chunk_rows.append(pd.read_csv(p.with_name(p.name.replace("_perturbations.txt", "_cells.tsv.gz")), sep="\t"))
    cohort = pd.concat(chunk_rows, ignore_index=True)
    control_label = json.loads((base / "chunk_runs/chunk_0000/run_meta.json").read_text())["control_label"]
    cohort = cohort[cohort.gene_target.isin(requested | {control_label})].drop_duplicates("cell_barcode")
    selected = pd.read_csv(base / "merged/selected_perturbed_cells_merged.tsv.gz", sep="\t")
    selected_barcodes = set(selected[selected.perturbation.isin(requested)].cell_barcode)
    raw_path = Path(old_meta["h5ad"]).resolve()
    with h5py.File(raw_path, "r") as h5:
        obs, var = read_elem(h5["obs"]), read_elem(h5["var"])
        rows = np.flatnonzero(obs.index.astype(str).isin(cohort.cell_barcode))
        assert len(rows) == len(cohort)
        labels = cohort.set_index("cell_barcode").gene_target.reindex(obs.index[rows]).astype(str)
        genes = np.arange(len(var))
        selector = "all measured features"
        x = h5["X"]
        if isinstance(x, h5py.Group) and x.attrs.get("encoding-type") == "csc_matrix":
            column = "n_cells_by_counts"
            top = np.argsort(var[column].to_numpy(), kind="stable")[-args.max_csc_features:]
            genes = np.union1d(top, np.flatnonzero(var.index.isin(targets)))
            selector = f"top {args.max_csc_features} atlas detection-count features plus requested targets"
        print(f"[cohort] {args.dataset}: {len(rows)} cells, {len(genes)} features", flush=True)
        raw = read_rows(h5, rows, genes)
        totals = obs.total_counts.iloc[rows].to_numpy() if "total_counts" in obs else np.asarray(raw.sum(axis=1)).ravel()
        if not np.isfinite(raw.data).all() or (raw.data < 0).any() or not np.allclose(raw.data, np.round(raw.data)):
            raise ValueError("Expected finite nonnegative integer raw counts")
        if (totals <= 0).any():
            raise ValueError("Zero library totals in cohort")
        normalized = sparse.diags(10000 / totals) @ raw
        normalized = normalized.tocsr().astype(np.float32)
        normalized.data = np.log1p(normalized.data)
        cohort_obs = obs.iloc[rows].copy()
        cohort_obs["target"] = labels.to_numpy()
        cohort_obs["historical_mixscape_selected"] = cohort_obs.index.isin(selected_barcodes)
        adata = ad.AnnData(X=normalized, obs=cohort_obs, var=var.iloc[genes].copy())
    sc.pp.highly_variable_genes(adata, n_top_genes=min(2000, adata.n_vars))
    # Prevent similarity from being defined primarily by each directly targeted gene.
    excluded = adata.var_names.isin(requested)
    mask = adata.var.highly_variable.to_numpy() & ~excluded
    adata.var["pca_features"] = mask
    sc.pp.pca(adata, n_comps=60, mask_var="pca_features", random_state=0)
    adata.X = adata.X.toarray()
    ms = pt.tl.Mixscape()
    print("[signature] shared dataset PCA and 30 nearest controls", flush=True)
    ms.perturbation_signature(
        adata, pert_key="target", control=control_label, n_neighbors=30,
        use_rep="X_pca", n_dims=60, batch_size=128, random_state=0,
    )
    # pertpy 1.1.1 stores reference minus cell. Orient as cell minus reference.
    residuals = -np.asarray(adata.layers["X_pert"])
    control_mask = adata.obs.target.eq(control_label).to_numpy()
    control_center = residuals[control_mask].mean(axis=0)
    residuals -= control_center
    features = ~excluded & (np.asarray(normalized.getnnz(axis=0)).ravel() >= 30)
    names = adata.var_names[features]
    assigned_profiles, selected_profiles, qc = [], [], []
    present_names, selected_names = [], []
    target_labels = adata.obs.target.to_numpy()
    ko = adata.obs.historical_mixscape_selected.to_numpy()
    for target in targets + randoms:
        match = target_labels == target
        selected_match = match & ko
        n_all, n_ko = int(match.sum()), int(selected_match.sum())
        half_agreement = np.nan
        if n_ko >= args.min_selected_cells:
            split = rng.permutation(np.flatnonzero(selected_match))
            first = residuals[split[:len(split)//2]][:, features].mean(axis=0)
            second = residuals[split[len(split)//2:]][:, features].mean(axis=0)
            if first.std() > 1e-8 and second.std() > 1e-8:
                half_agreement = float(np.corrcoef(first, second)[0, 1])
        qc.append({
            "perturbation": target, "category": "requested" if target in targets else "random_comparator",
            "is_goi": target in goi, "n_assigned_cells": n_all, "n_selected_cells": n_ko,
            "n_profile_cells": n_ko,
            "split_half_profile_correlation": half_agreement,
            "status": (
                "clustered" if n_ko >= args.min_selected_cells else
                "below_min_selected_cells" if n_ko else
                "no_selected_cells" if n_all else "not_in_processed_cohort"
            ),
        })
        if n_all >= args.min_selected_cells:
            present_names.append(target)
            assigned_profiles.append(residuals[match][:, features].mean(axis=0))
        if n_ko >= args.min_selected_cells:
            selected_names.append(target)
            selected_profiles.append(residuals[selected_match][:, features].mean(axis=0))
    qc = pd.DataFrame(qc)
    qc.to_csv(args.outdir / "target_coverage_qc.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": adata.var_names, "clustering_feature": features, "direct_target_excluded": excluded}).to_csv(
        args.outdir / "feature_manifest.tsv", sep="\t", index=False
    )
    cohort_obs[["target", "historical_mixscape_selected"]].to_csv(args.outdir / "cohort_cells.tsv.gz", sep="\t")
    summaries = {}
    for cohort_name, vectors, row_names in [
        ("selected_cells", selected_profiles, selected_names),
        ("all_assigned_cells_exploratory", assigned_profiles, present_names),
    ]:
        out = args.outdir / cohort_name
        out.mkdir(exist_ok=True)
        frame = pd.DataFrame(vectors, index=row_names, columns=names)
        frame.index.name = "perturbation"
        cohort_qc = qc.copy()
        if cohort_name.startswith("all"):
            cohort_qc["n_profile_cells"] = cohort_qc.n_assigned_cells
        summaries[cohort_name] = save_clustering(frame, cohort_qc, out, f"{args.dataset}: {cohort_name.replace('_', ' ')}")
    # Comparable raw-expression correlations on precisely the same rows/features.
    order = selected_names
    means = np.array([adata.X[(target_labels == t) & ko][:, features].mean(axis=0) for t in order])
    raw_corr = np.corrcoef(means) if len(order) > 1 else None
    if raw_corr is not None:
        pd.DataFrame(raw_corr, index=order, columns=order).to_csv(args.outdir / "matched_absolute_expression_correlations.tsv", sep="\t")
        summaries["matched_absolute_expression_median_correlation"] = float(np.median(raw_corr[np.triu_indices(len(order), 1)]))
    meta = {
        "dataset": args.dataset, "method": "recomputed pertpy nearest-control X_pert, historical KO labels retained",
        "is_exact_historical_signature_reproduction": False,
        "reference": "one shared PCA/control cohort per dataset; controls from original manifests",
        "feature_selector": selector, "n_input_features": adata.n_vars, "n_clustering_features": int(features.sum()),
        "direct_targets_excluded_from_clustering": sorted(requested),
        "normalization": "10000/full-cell raw library total; log1p",
        "sign": "cell minus nearest-control reference, centered by mean control residual",
        "n_controls": int(control_mask.sum()), "n_cells": adata.n_obs, "min_profile_cells": args.min_selected_cells,
        "random_seed": 0, "n_neighbors": 30, "pca_dimensions": 60, "n_random_comparators": len(randoms),
        "raw_h5ad": str(raw_path), "raw_size": raw_path.stat().st_size, "raw_mtime_ns": raw_path.stat().st_mtime_ns,
        "target_list_sha256": hashlib.sha256(target_file.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "pertpy_version": pt.__version__, "scanpy_version": sc.__version__,
        "elapsed_seconds": time.time() - started, "summaries": summaries,
        "interpretation": "Exploratory descriptive clusters. Historical cell selection and small cohorts can bias or destabilize profiles. No biological-replicate inference.",
    }
    (args.outdir / "run_meta.json").write_text(json.dumps(meta, indent=2))
    print("[done]", json.dumps(summaries), flush=True)


if __name__ == "__main__":
    main()
