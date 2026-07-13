#!/usr/bin/env python3
"""Target-vs-random clustering with two arms:
1) pseudobulk expression from KO-assigned cells (Mixscape + Mixscale),
2) DEG profile from all configured DEG streams (non-significant values set to 0).
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def _read_target_list(path: Path) -> list[str]:
    rows: list[str] = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        rows.append(line)
    return list(dict.fromkeys(rows))


def _to_dense_row(x) -> np.ndarray:
    if sparse.issparse(x):
        return np.asarray(x).ravel()
    return np.asarray(x).ravel()


def _cluster_matrix(
    matrix: pd.DataFrame,
    n_clusters: int,
    linkage_method: str,
    distance_metric: str,
) -> tuple[pd.DataFrame, list[str]]:
    if matrix.empty:
        return pd.DataFrame(columns=["perturbation", "cluster", "order_index"]), []

    perts = matrix.index.astype(str).tolist()
    if matrix.shape[0] <= 1 or matrix.shape[1] == 0:
        out = pd.DataFrame(
            {
                "perturbation": perts,
                "cluster": np.ones(len(perts), dtype=int),
                "order_index": np.arange(len(perts), dtype=int),
            }
        )
        return out, perts

    vals = matrix.to_numpy(dtype=float)
    mu = vals.mean(axis=0, keepdims=True)
    sd = vals.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    vals = (vals - mu) / sd

    d = pdist(vals, metric=str(distance_metric))
    d = np.nan_to_num(d, nan=0.0, posinf=0.0, neginf=0.0)
    z = linkage(d, method=str(linkage_method))
    leaves = leaves_list(z).astype(int).tolist()
    ordered = [perts[i] for i in leaves]
    labels = fcluster(z, t=max(1, min(int(n_clusters), len(perts))), criterion="maxclust").astype(int)

    order_map = {p: i for i, p in enumerate(ordered)}
    out = pd.DataFrame({"perturbation": perts, "cluster": labels})
    out["order_index"] = out["perturbation"].map(order_map).astype(int)
    out = out.sort_values("order_index").reset_index(drop=True)
    return out, ordered


def _pairwise_similarity(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix.empty:
        return pd.DataFrame(index=matrix.index, columns=matrix.index)
    vals = matrix.to_numpy(dtype=float)
    if vals.shape[0] <= 1:
        sim = np.ones((vals.shape[0], vals.shape[0]), dtype=float)
    else:
        mu = vals.mean(axis=0, keepdims=True)
        sd = vals.std(axis=0, keepdims=True)
        sd[sd == 0] = 1.0
        vals = (vals - mu) / sd
        d = squareform(pdist(vals, metric="euclidean"))
        sim = 1.0 / (1.0 + d)
    out = pd.DataFrame(sim, index=matrix.index.astype(str), columns=matrix.index.astype(str))
    out.index.name = "perturbation"
    return out


def _select_heatmap_columns(matrix: pd.DataFrame, max_genes: int, gene_order: str) -> pd.DataFrame:
    if matrix.empty:
        return matrix
    out = matrix
    if int(max_genes) > 0 and out.shape[1] > int(max_genes):
        if str(gene_order) == "absolute_mean":
            score = out.abs().mean(axis=0)
        else:
            score = out.var(axis=0)
        keep = score.sort_values(ascending=False).head(int(max_genes)).index.tolist()
        out = out[keep]
    return out


def _group_color(group: str) -> str:
    if group == "target":
        return "#b2182b"
    if group == "random_control":
        return "#2166ac"
    if group == "control":
        return "#1b7837"
    return "#000000"


def _write_pairwise_heatmap(
    matrix: pd.DataFrame,
    order: list[str],
    groups: dict[str, str],
    out_png: Path,
    title: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    if matrix.empty or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_axis_off()
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(out_png, dpi=180)
        plt.close(fig)
        return

    ordered = [p for p in order if p in matrix.index and p in matrix.columns]
    if not ordered:
        ordered = matrix.index.astype(str).tolist()
    mat = matrix.loc[ordered, ordered]

    n = mat.shape[0]
    fig_size = min(20.0, max(6.0, 0.28 * n + 2.0))
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="equal", cmap="viridis", vmin=0.0, vmax=1.0, interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("Perturbations")
    ax.set_ylabel("Perturbations")

    if n <= 120:
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))
        ax.set_xticklabels(ordered, rotation=90, fontsize=7)
        ax.set_yticklabels(ordered, fontsize=7)
        for tick in ax.get_yticklabels():
            tick.set_color(_group_color(groups.get(tick.get_text(), "")))
        for tick in ax.get_xticklabels():
            tick.set_color(_group_color(groups.get(tick.get_text(), "")))
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Similarity")
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def _write_volcano_plot(
    df: pd.DataFrame,
    perturbation: str,
    recurrent_genes: set[str],
    min_abs_effect: float,
    max_pval_adj: float,
    out_png: Path,
) -> dict[str, object]:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.text(0.5, 0.5, f"No volcano data for {perturbation}", ha="center", va="center")
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(out_png, dpi=180)
        plt.close(fig)
        return {"n_points": 0, "n_sig": 0, "n_recurrent_hits": 0}

    x = pd.to_numeric(df["effect_value"], errors="coerce")
    p = pd.to_numeric(df["pval_adj"], errors="coerce")
    valid = x.notna() & p.notna()
    d = df.loc[valid, ["gene"]].copy()
    d["x"] = x.loc[valid].to_numpy(dtype=float)
    d["p"] = p.loc[valid].clip(lower=1e-300).to_numpy(dtype=float)
    d["y"] = -np.log10(d["p"])
    d["sig"] = (np.abs(d["x"]) >= float(min_abs_effect)) & (d["p"] <= float(max_pval_adj))
    d["recurrent"] = d["gene"].astype(str).isin(recurrent_genes)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(d.loc[~d["sig"], "x"], d.loc[~d["sig"], "y"], s=8, c="#bdbdbd", alpha=0.6, linewidths=0)
    ax.scatter(d.loc[d["sig"], "x"], d.loc[d["sig"], "y"], s=10, c="#ef3b2c", alpha=0.8, linewidths=0)
    ax.scatter(d.loc[d["recurrent"], "x"], d.loc[d["recurrent"], "y"], s=18, c="#08519c", alpha=0.9, linewidths=0)
    ax.axvline(float(min_abs_effect), color="#666666", linestyle="--", linewidth=0.8)
    ax.axvline(-float(min_abs_effect), color="#666666", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(max(float(max_pval_adj), 1e-300)), color="#666666", linestyle="--", linewidth=0.8)
    ax.set_title(f"{perturbation}: Volcano")
    ax.set_xlabel("Effect")
    ax.set_ylabel("-log10(adj p)")

    # Label recurrent genes and the perturbation gene itself if present.
    label_df = d[d["recurrent"]].copy()
    if perturbation in set(d["gene"].astype(str).tolist()):
        label_df = pd.concat([label_df, d[d["gene"].astype(str) == perturbation]], axis=0, ignore_index=True).drop_duplicates("gene")
    label_df = label_df.sort_values(["y", "x"], ascending=[False, False]).head(25)
    for _, r in label_df.iterrows():
        ax.text(float(r["x"]), float(r["y"]), str(r["gene"]), fontsize=7)

    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)

    return {
        "n_points": int(d.shape[0]),
        "n_sig": int(d["sig"].sum()),
        "n_recurrent_hits": int(d["recurrent"].sum()),
    }


def _standardize_de_stream(path: Path, stream: str) -> tuple[pd.DataFrame, dict[str, object]]:
    preferred_path = path
    if path.name == "perturbation_differential_genes.tsv.gz":
        full_path = path.with_name("perturbation_differential_genes_full.tsv.gz")
        if full_path.exists():
            preferred_path = full_path

    meta = {"stream": stream, "path": str(preferred_path), "loaded_rows": 0, "used_rows": 0, "effect_mode": "none"}
    if not preferred_path.exists():
        cols = ["perturbation", "gene", "effect_value", "pval_adj", "rank", "stream", "effect_mode"]
        return pd.DataFrame(columns=cols), meta

    df = pd.read_csv(preferred_path, sep="\t", compression="infer")
    meta["loaded_rows"] = int(df.shape[0])
    if df.empty or "perturbation" not in df.columns or "gene" not in df.columns:
        cols = ["perturbation", "gene", "effect_value", "pval_adj", "rank", "stream", "effect_mode"]
        return pd.DataFrame(columns=cols), meta

    out = df.copy()
    out["perturbation"] = out["perturbation"].astype(str)
    out["gene"] = out["gene"].astype(str)
    out["rank"] = pd.to_numeric(out.get("rank", np.nan), errors="coerce")
    out["pval_adj"] = pd.to_numeric(out.get("pval_adj", out.get("pval", np.nan)), errors="coerce")

    effect_mode = "none"
    effect = None
    if "logfoldchange" in out.columns:
        v = pd.to_numeric(out["logfoldchange"], errors="coerce")
        if bool(v.notna().any()):
            effect = v
            effect_mode = "logfoldchange"
    if effect is None and "score" in out.columns:
        v = pd.to_numeric(out["score"], errors="coerce")
        if bool(v.notna().any()):
            effect = v
            effect_mode = "score"
    if effect is None and "p_weight" in out.columns:
        v = pd.to_numeric(out["p_weight"], errors="coerce")
        if bool(v.notna().any()):
            effect = v
            effect_mode = "p_weight"
    if effect is None:
        v = pd.to_numeric(out.get("rank", np.nan), errors="coerce")
        if bool(v.notna().any()):
            effect = 1.0 / v.clip(lower=1.0)
            effect_mode = "rank_inverse"

    if effect is None:
        cols = ["perturbation", "gene", "effect_value", "pval_adj", "rank", "stream", "effect_mode"]
        return pd.DataFrame(columns=cols), meta

    out["effect_value"] = pd.to_numeric(effect, errors="coerce")
    out = out.dropna(subset=["effect_value"]).copy()
    out["stream"] = str(stream)
    out["effect_mode"] = effect_mode
    meta["effect_mode"] = effect_mode
    meta["used_rows"] = int(out.shape[0])
    return out[["perturbation", "gene", "effect_value", "pval_adj", "rank", "stream", "effect_mode"]], meta


def _combine_deg_streams(
    de_tables: list[Path],
    stream_names: list[str],
    min_abs_effect: float,
    max_pval_adj: float,
    rank_significant_max: int,
) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, object]], int]:
    def _aggregate_standardized(raw_df: pd.DataFrame, denom: int) -> pd.DataFrame:
        if raw_df.empty:
            cols = ["perturbation", "gene", "effect_for_cluster", "abs_effect", "n_streams_detected"]
            return pd.DataFrame(columns=cols)
        agg_df = (
            raw_df.groupby(["perturbation", "gene"], as_index=False)
            .agg(
                sum_effect_for_cluster=("effect_for_cluster", "sum"),
                abs_effect=("abs_effect", "max"),
                n_streams_detected=("stream", "nunique"),
            )
        )
        denom_val = max(1, int(denom))
        agg_df["effect_for_cluster"] = agg_df["sum_effect_for_cluster"] / float(denom_val)
        return agg_df.drop(columns=["sum_effect_for_cluster"])

    frames: list[pd.DataFrame] = []
    metas: list[dict[str, object]] = []
    for path, stream in zip(de_tables, stream_names):
        df, meta = _standardize_de_stream(path, stream)
        frames.append(df)
        metas.append(meta)

    raw = pd.concat(frames, axis=0, ignore_index=True) if frames else pd.DataFrame()
    n_streams_total = int(len(stream_names))
    if raw.empty:
        cols = ["perturbation", "gene", "effect_for_cluster", "abs_effect", "n_streams_detected"]
        return pd.DataFrame(columns=cols), raw, metas, n_streams_total

    sig = np.ones(raw.shape[0], dtype=bool)
    has_p = raw["pval_adj"].notna().to_numpy()
    if float(max_pval_adj) >= 0:
        sig[has_p] &= raw.loc[has_p, "pval_adj"].to_numpy(dtype=float) <= float(max_pval_adj)

    is_rank_mode = raw["effect_mode"].eq("rank_inverse").to_numpy()
    if int(rank_significant_max) > 0:
        rank_ok = raw["rank"].notna().to_numpy()
        sig[is_rank_mode & rank_ok] &= (
            raw.loc[is_rank_mode & rank_ok, "rank"].to_numpy(dtype=float) <= float(rank_significant_max)
        )
        sig[is_rank_mode & ~rank_ok] = False

    non_rank = ~is_rank_mode
    if float(min_abs_effect) > 0:
        sig[non_rank] &= np.abs(raw.loc[non_rank, "effect_value"].to_numpy(dtype=float)) >= float(min_abs_effect)

    raw["effect_for_cluster"] = np.where(sig, raw["effect_value"].to_numpy(dtype=float), 0.0)
    raw["abs_effect"] = np.abs(raw["effect_value"].to_numpy(dtype=float))
    agg = _aggregate_standardized(raw, denom=n_streams_total)
    return agg, raw, metas, n_streams_total


def _aggregate_standardized_deg(raw_df: pd.DataFrame, denom: int) -> pd.DataFrame:
    if raw_df.empty:
        cols = ["perturbation", "gene", "effect_for_cluster", "abs_effect", "n_streams_detected"]
        return pd.DataFrame(columns=cols)
    agg_df = (
        raw_df.groupby(["perturbation", "gene"], as_index=False)
        .agg(
            sum_effect_for_cluster=("effect_for_cluster", "sum"),
            abs_effect=("abs_effect", "max"),
            n_streams_detected=("stream", "nunique"),
        )
    )
    denom_val = max(1, int(denom))
    agg_df["effect_for_cluster"] = agg_df["sum_effect_for_cluster"] / float(denom_val)
    return agg_df.drop(columns=["sum_effect_for_cluster"])


def _seed_for_label(base_seed: int, label: str) -> int:
    offset = sum((i + 1) * ord(ch) for i, ch in enumerate(label))
    return int((int(base_seed) + offset) % (2**32 - 1))


def _sample_target_and_random(
    requested_targets: list[str],
    candidate_perts: set[str],
    random_match_target_count: bool,
    n_random_controls: int,
    rng: np.random.Generator,
) -> tuple[list[str], list[str], list[str]]:
    targets_present = [t for t in requested_targets if t in candidate_perts]
    missing_targets = [t for t in requested_targets if t not in set(targets_present)]
    candidate_pool = sorted(set(candidate_perts) - set(targets_present))
    if bool(random_match_target_count):
        desired_random = len(requested_targets)
    else:
        desired_random = int(n_random_controls) if int(n_random_controls) > 0 else len(targets_present)
    n_random = min(int(desired_random), len(candidate_pool))
    sampled_random = (
        sorted(rng.choice(np.asarray(candidate_pool, dtype=object), size=n_random, replace=False).tolist())
        if n_random > 0
        else []
    )
    return targets_present, missing_targets, sampled_random


def _build_selected_df(
    selected_ko_perts: list[str],
    targets_present: list[str],
    sampled_random: list[str],
    control_label: str,
    present_in_expression: set[str] | None = None,
    present_in_de: set[str] | None = None,
    used_counts: dict[str, int] | None = None,
) -> pd.DataFrame:
    if present_in_expression is None:
        present_in_expression = set()
    if present_in_de is None:
        present_in_de = set()
    if used_counts is None:
        used_counts = {}
    control = str(control_label)
    return pd.DataFrame(
        {
            "perturbation": selected_ko_perts + [control],
            "group": (["target"] * len(targets_present)) + (["random_control"] * len(sampled_random)) + ["control"],
            "is_target": ([True] * len(targets_present)) + ([False] * len(sampled_random)) + [False],
            "is_control": ([False] * len(selected_ko_perts)) + [True],
            "present_in_de": ([p in present_in_de for p in selected_ko_perts] + [True]),
            "present_in_expression": ([p in present_in_expression for p in selected_ko_perts] + [control in present_in_expression]),
            "n_ko_cells_used": ([int(used_counts.get(p, 0)) for p in selected_ko_perts] + [int(used_counts.get(control, 0))]),
        }
    )


def _load_mixscape_selected(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["cell_barcode", "perturbation"])
    df = pd.read_csv(path, sep="\t", compression="infer", dtype="string")
    if df.empty or "cell_barcode" not in df.columns or "perturbation" not in df.columns:
        return pd.DataFrame(columns=["cell_barcode", "perturbation"])
    out = df[["cell_barcode", "perturbation"]].dropna().copy()
    out["cell_barcode"] = out["cell_barcode"].astype(str)
    out["perturbation"] = out["perturbation"].astype(str)
    return out.drop_duplicates()


def _collect_mixscale_assigned(path: Path, score_threshold: float, selected_perts: set[str] | None = None) -> pd.DataFrame:
    cols = ["cell_barcode", "perturbation", "mixscale_score"]
    if not path.exists():
        return pd.DataFrame(columns=["cell_barcode", "perturbation"])
    out_chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, sep="\t", compression="infer", usecols=cols, chunksize=1_000_000, dtype="string"):
        c = chunk.copy()
        c["cell_barcode"] = c["cell_barcode"].astype(str)
        c["perturbation"] = c["perturbation"].astype(str)
        c["mixscale_score"] = pd.to_numeric(c["mixscale_score"], errors="coerce")
        c = c[c["mixscale_score"] > float(score_threshold)]
        if selected_perts is not None:
            c = c[c["perturbation"].isin(selected_perts)]
        if c.empty:
            continue
        out_chunks.append(c[["cell_barcode", "perturbation"]])
    if not out_chunks:
        return pd.DataFrame(columns=["cell_barcode", "perturbation"])
    out = pd.concat(out_chunks, axis=0, ignore_index=True)
    return out.drop_duplicates()


def _build_assignment_map(
    selected_perts: list[str],
    mixscape_df: pd.DataFrame,
    mixscale_df: pd.DataFrame,
    logic: str,
) -> dict[str, set[str]]:
    mixscape_map = {p: set(g["cell_barcode"].astype(str).tolist()) for p, g in mixscape_df.groupby("perturbation", sort=False)}
    mixscale_map = {p: set(g["cell_barcode"].astype(str).tolist()) for p, g in mixscale_df.groupby("perturbation", sort=False)}
    out: dict[str, set[str]] = {}
    for p in selected_perts:
        a = mixscape_map.get(p, set())
        b = mixscale_map.get(p, set())
        cells = (a & b) if logic == "intersection" else (a | b)
        if cells:
            out[p] = cells
    return out


def _load_control_barcodes(
    h5ad: Path,
    pert_col: str,
    control_label: str,
    max_control_cells: int,
    seed: int,
) -> list[str]:
    adata = sc.read_h5ad(h5ad, backed="r")
    obs = adata.obs[[pert_col]].copy()
    obs.index = obs.index.astype(str)
    obs[pert_col] = obs[pert_col].astype(str)
    try:
        ctrl = obs.index[obs[pert_col] == str(control_label)].astype(str).tolist()
    finally:
        adata.file.close()
    if int(max_control_cells) > 0 and len(ctrl) > int(max_control_cells):
        rng = np.random.default_rng(int(seed))
        ctrl = rng.choice(np.asarray(ctrl, dtype=object), size=int(max_control_cells), replace=False).tolist()
    return list(ctrl)


def _build_pseudobulk_expression(
    h5ad: Path,
    assignment_map: dict[str, set[str]],
    control_label: str,
    control_barcodes: list[str],
    normalize_target_sum: float,
    log1p: bool,
) -> tuple[pd.DataFrame, dict[str, int]]:
    groups = list(assignment_map.keys())
    keep = set(control_barcodes)
    for bars in assignment_map.values():
        keep.update(bars)
    if not keep:
        return pd.DataFrame(index=pd.Index(groups + [control_label], name="perturbation")), {}

    adata_backed = sc.read_h5ad(h5ad, backed="r")
    obs_names = pd.Index(adata_backed.obs_names.astype(str))
    keep_mask = obs_names.isin(keep)
    adata = adata_backed[keep_mask].to_memory()
    adata_backed.file.close()

    adata.obs_names = adata.obs_names.astype(str)
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    genes = adata.var_names.astype(str).tolist()

    rows = []
    row_names = []
    used_counts: dict[str, int] = {}
    obs_index = pd.Index(adata.obs_names.astype(str))

    for p in groups:
        bars = list(assignment_map.get(p, set()))
        if not bars:
            continue
        mask = obs_index.isin(bars)
        n = int(mask.sum())
        if n == 0:
            continue
        rows.append(_to_dense_row(X[mask].sum(axis=0)).astype(float))
        row_names.append(p)
        used_counts[p] = n

    ctrl_mask = obs_index.isin(control_barcodes)
    n_ctrl = int(ctrl_mask.sum())
    if n_ctrl > 0:
        rows.append(_to_dense_row(X[ctrl_mask].sum(axis=0)).astype(float))
        row_names.append(str(control_label))
        used_counts[str(control_label)] = n_ctrl

    if not rows:
        return pd.DataFrame(index=pd.Index(groups + [control_label], name="perturbation")), used_counts

    mat = np.vstack(rows).astype(float)
    if float(normalize_target_sum) > 0:
        rs = mat.sum(axis=1, keepdims=True)
        rs[rs <= 0] = 1.0
        mat = mat * (float(normalize_target_sum) / rs)
    if bool(log1p):
        mat = np.log1p(mat)

    df = pd.DataFrame(mat, index=row_names, columns=genes)
    df.index.name = "perturbation"
    return df, used_counts


def _build_deg_matrix(
    de_agg: pd.DataFrame,
    selected_perts: list[str],
    control_label: str,
    top_genes_per_perturbation: int,
    min_feature_genes: int,
) -> pd.DataFrame:
    if de_agg.empty or not selected_perts:
        return pd.DataFrame(index=pd.Index(selected_perts + [str(control_label)], name="perturbation"))

    sub = de_agg[de_agg["perturbation"].isin(set(selected_perts))].copy()
    if sub.empty:
        return pd.DataFrame(index=pd.Index(selected_perts + [str(control_label)], name="perturbation"))

    feature_genes: set[str] = set()
    if int(top_genes_per_perturbation) > 0:
        for p in selected_perts:
            g = sub.loc[sub["perturbation"] == p, ["gene", "abs_effect"]]
            if g.empty:
                continue
            feature_genes.update(g.nlargest(int(top_genes_per_perturbation), "abs_effect")["gene"].astype(str).tolist())
    else:
        feature_genes.update(sub["gene"].astype(str).tolist())

    if int(min_feature_genes) > 0 and len(feature_genes) < int(min_feature_genes):
        global_top = (
            sub.groupby("gene", as_index=False)["abs_effect"]
            .max()
            .sort_values("abs_effect", ascending=False)
            .head(int(min_feature_genes))["gene"]
            .astype(str)
            .tolist()
        )
        feature_genes.update(global_top)

    genes = sorted(feature_genes)
    if not genes:
        return pd.DataFrame(index=pd.Index(selected_perts + [str(control_label)], name="perturbation"))

    sub = sub[sub["gene"].isin(genes)].copy()
    mat = sub.pivot_table(
        index="perturbation",
        columns="gene",
        values="effect_for_cluster",
        aggfunc="mean",
        fill_value=0.0,
    )
    mat = mat.reindex(selected_perts).fillna(0.0)
    mat.loc[str(control_label)] = 0.0
    mat = mat.reindex(selected_perts + [str(control_label)]).fillna(0.0)
    mat.index.name = "perturbation"
    return mat


def _subset_ora_for_targets(
    ora_tables: list[Path],
    ora_streams: list[str],
    targets: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for path, stream in zip(ora_tables, ora_streams):
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t", compression="infer")
        if df.empty or "perturbation" not in df.columns:
            continue
        sub = df[df["perturbation"].astype(str).isin(set(targets))].copy()
        if sub.empty:
            continue
        sub.insert(0, "stream", str(stream))
        rows.append(sub)
    if not rows:
        empty_terms = pd.DataFrame(columns=["stream", "perturbation", "direction", "term", "p_adj_bh", "p_value"])
        empty_summary = pd.DataFrame(columns=["stream", "perturbation", "direction", "n_terms", "best_term", "best_p_adj_or_p"])
        return empty_terms, empty_summary

    terms = pd.concat(rows, axis=0, ignore_index=True)
    pcol = "p_adj_bh" if "p_adj_bh" in terms.columns else ("p_adj" if "p_adj" in terms.columns else ("p_value" if "p_value" in terms.columns else None))
    if pcol is None:
        terms["_best_p"] = np.nan
    else:
        terms["_best_p"] = pd.to_numeric(terms[pcol], errors="coerce")
    grp_cols = [c for c in ["stream", "perturbation", "direction"] if c in terms.columns]
    summary_rows = []
    for keys, g in terms.groupby(grp_cols, sort=True):
        g2 = g.sort_values("_best_p", ascending=True, na_position="last")
        key_vals = keys if isinstance(keys, tuple) else (keys,)
        rec = dict(zip(grp_cols, key_vals))
        rec["n_terms"] = int(g2.shape[0])
        rec["best_term"] = str(g2.iloc[0]["term"]) if ("term" in g2.columns and g2.shape[0] > 0) else np.nan
        rec["best_p_adj_or_p"] = float(g2.iloc[0]["_best_p"]) if g2.shape[0] > 0 and pd.notna(g2.iloc[0]["_best_p"]) else np.nan
        summary_rows.append(rec)
    summary = pd.DataFrame(summary_rows)
    return terms.drop(columns=["_best_p"]), summary


def run_analysis(args: argparse.Namespace) -> None:
    args.de_tables = [Path(str(p)) for p in args.de_tables]
    args.ora_tables = [Path(str(p)) for p in args.ora_tables]
    args.mixscape_selected = Path(args.mixscape_selected)
    args.mixscale_cells = Path(args.mixscale_cells)
    args.chunk_summary = Path(args.chunk_summary)
    args.h5ad = Path(args.h5ad)
    args.target_list = Path(args.target_list)
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    requested_targets = _read_target_list(args.target_list)
    if not requested_targets:
        raise RuntimeError(f"Target list is empty: {args.target_list}")

    stream_names = [s.strip() for s in str(args.deg_streams).split(",") if s.strip()]
    if len(stream_names) != len(args.de_tables):
        stream_names = [f"stream_{i+1}" for i in range(len(args.de_tables))]
    ora_streams = [s.strip() for s in str(args.ora_streams).split(",") if s.strip()]
    if len(ora_streams) != len(args.ora_tables):
        ora_streams = [f"stream_{i+1}" for i in range(len(args.ora_tables))]

    de_agg, de_raw_all, de_stream_meta, n_streams_total = _combine_deg_streams(
        de_tables=args.de_tables,
        stream_names=stream_names,
        min_abs_effect=float(args.min_abs_effect),
        max_pval_adj=float(args.max_pval_adj),
        rank_significant_max=int(args.rank_significant_max),
    )
    de_perts = set(de_agg["perturbation"].astype(str).tolist()) if not de_agg.empty else set()
    de_agg_by_stream: dict[str, pd.DataFrame] = {}
    if not de_raw_all.empty and "stream" in de_raw_all.columns:
        for s in stream_names:
            raw_s = de_raw_all[de_raw_all["stream"].astype(str) == str(s)].copy()
            de_agg_by_stream[str(s)] = _aggregate_standardized_deg(raw_s, denom=1)
    else:
        for s in stream_names:
            de_agg_by_stream[str(s)] = _aggregate_standardized_deg(pd.DataFrame(), denom=1)

    mixscape_all = _load_mixscape_selected(args.mixscape_selected)
    mixscale_all_assigned = _collect_mixscale_assigned(
        path=args.mixscale_cells,
        score_threshold=float(args.mixscale_score_threshold),
        selected_perts=None,
    )
    mixscape_map = {
        p: set(g["cell_barcode"].astype(str).tolist()) for p, g in mixscape_all.groupby("perturbation", sort=False)
    } if not mixscape_all.empty else {}
    mixscale_map = {
        p: set(g["cell_barcode"].astype(str).tolist()) for p, g in mixscale_all_assigned.groupby("perturbation", sort=False)
    } if not mixscale_all_assigned.empty else {}
    mixscape_perts = set(mixscape_all["perturbation"].astype(str).tolist())
    mixscale_perts = set(mixscale_all_assigned["perturbation"].astype(str).tolist())
    expr_candidate_perts = (mixscape_perts & mixscale_perts) if str(args.assignment_logic) == "intersection" else (mixscape_perts | mixscale_perts)

    rng = np.random.default_rng(_seed_for_label(int(args.random_seed), "legacy_combined"))
    targets_present, missing_targets, sampled_random = _sample_target_and_random(
        requested_targets=requested_targets,
        candidate_perts=(de_perts & expr_candidate_perts),
        random_match_target_count=bool(args.random_match_target_count),
        n_random_controls=int(args.n_random_controls),
        rng=rng,
    )

    selected_ko_perts = list(targets_present) + list(sampled_random)
    mixscape_sel = mixscape_all[mixscape_all["perturbation"].isin(set(selected_ko_perts))].copy()
    mixscale_sel = _collect_mixscale_assigned(
        path=args.mixscale_cells,
        score_threshold=float(args.mixscale_score_threshold),
        selected_perts=set(selected_ko_perts),
    )
    assignment_map = _build_assignment_map(
        selected_perts=selected_ko_perts,
        mixscape_df=mixscape_sel,
        mixscale_df=mixscale_sel,
        logic=str(args.assignment_logic),
    )
    selected_ko_perts = [p for p in selected_ko_perts if p in assignment_map]

    control_barcodes = _load_control_barcodes(
        h5ad=args.h5ad,
        pert_col=args.pert_col,
        control_label=args.control_label,
        max_control_cells=int(args.max_control_cells),
        seed=int(args.random_seed),
    )

    expr_mat, expr_used_counts = _build_pseudobulk_expression(
        h5ad=args.h5ad,
        assignment_map=assignment_map,
        control_label=str(args.control_label),
        control_barcodes=control_barcodes,
        normalize_target_sum=float(args.expression_normalize_target_sum),
        log1p=bool(args.expression_log1p),
    )

    deg_mat = _build_deg_matrix(
        de_agg=de_agg,
        selected_perts=selected_ko_perts,
        control_label=str(args.control_label),
        top_genes_per_perturbation=int(args.top_genes_per_perturbation),
        min_feature_genes=int(args.min_feature_genes),
    )

    selected_df = _build_selected_df(
        selected_ko_perts=selected_ko_perts,
        targets_present=targets_present,
        sampled_random=sampled_random,
        control_label=str(args.control_label),
        present_in_expression=set(expr_mat.index.astype(str).tolist()),
        present_in_de=set(deg_mat.index.astype(str).tolist()),
        used_counts=expr_used_counts,
    )
    row_groups = dict(zip(selected_df["perturbation"].astype(str), selected_df["group"].astype(str)))

    expr_clusters, expr_ordered = _cluster_matrix(
        matrix=expr_mat,
        n_clusters=int(args.n_clusters),
        linkage_method=str(args.linkage),
        distance_metric=str(args.distance),
    )
    deg_clusters, deg_ordered = _cluster_matrix(
        matrix=deg_mat,
        n_clusters=int(args.n_clusters),
        linkage_method=str(args.linkage),
        distance_metric=str(args.distance),
    )
    expr_clusters = expr_clusters.merge(selected_df[["perturbation", "group", "is_target", "is_control"]], on="perturbation", how="left")
    deg_clusters = deg_clusters.merge(selected_df[["perturbation", "group", "is_target", "is_control"]], on="perturbation", how="left")

    expr_pair = _pairwise_similarity(expr_mat)
    deg_pair = _pairwise_similarity(deg_mat)

    expr_heatmap_mat = _select_heatmap_columns(
        matrix=expr_mat,
        max_genes=int(args.heatmap_max_genes),
        gene_order=str(args.heatmap_gene_order),
    )
    deg_heatmap_mat = _select_heatmap_columns(
        matrix=deg_mat,
        max_genes=int(args.heatmap_max_genes),
        gene_order=str(args.heatmap_gene_order),
    )

    # Volcano plots for targets: all if <N, else random N targets.
    volc_targets = list(targets_present)
    if len(volc_targets) > int(args.volcano_max_targets):
        volc_targets = sorted(
            rng.choice(np.asarray(volc_targets, dtype=object), size=int(args.volcano_max_targets), replace=False).tolist()
        )
    recurrent_cut = max(1, int(math.ceil(float(args.recurrent_gene_fraction) * max(1, len(targets_present)))))
    recurrent_rows = de_raw_all[
        de_raw_all["perturbation"].astype(str).isin(set(targets_present))
        & (pd.to_numeric(de_raw_all["effect_for_cluster"], errors="coerce").fillna(0.0) != 0.0)
    ][["perturbation", "gene"]].drop_duplicates()
    recurrent_counts = recurrent_rows.groupby("gene", as_index=False)["perturbation"].nunique() if not recurrent_rows.empty else pd.DataFrame(columns=["gene", "perturbation"])
    recurrent_genes = set(
        recurrent_counts.loc[recurrent_counts["perturbation"] >= recurrent_cut, "gene"].astype(str).tolist()
    ) if not recurrent_counts.empty else set()
    volc_stream_priority = [s.strip() for s in str(args.volcano_stream_priority).split(",") if s.strip()]

    volcano_dir = args.outdir / "volcano"
    volcano_rows = []
    for pert in volc_targets:
        chosen = pd.DataFrame()
        chosen_stream = ""
        for s in volc_stream_priority:
            cand = de_raw_all[
                (de_raw_all["perturbation"].astype(str) == str(pert))
                & (de_raw_all["stream"].astype(str) == str(s))
                & (de_raw_all["effect_mode"].astype(str).isin(["logfoldchange", "score"]))
                & pd.to_numeric(de_raw_all["pval_adj"], errors="coerce").notna()
            ].copy()
            if not cand.empty:
                chosen = cand
                chosen_stream = s
                break
        if chosen.empty:
            chosen = de_raw_all[
                (de_raw_all["perturbation"].astype(str) == str(pert))
                & (de_raw_all["effect_mode"].astype(str).isin(["logfoldchange", "score"]))
                & pd.to_numeric(de_raw_all["pval_adj"], errors="coerce").notna()
            ].copy()
            if not chosen.empty:
                chosen_stream = str(chosen["stream"].iloc[0])

        plot_path = volcano_dir / f"volcano_{str(pert).replace('/', '_')}.png"
        stats = _write_volcano_plot(
            df=chosen,
            perturbation=str(pert),
            recurrent_genes=recurrent_genes,
            min_abs_effect=float(args.min_abs_effect),
            max_pval_adj=float(args.max_pval_adj),
            out_png=plot_path,
        )
        volcano_rows.append(
            {
                "perturbation": str(pert),
                "stream_used": chosen_stream,
                "plot_path": str(plot_path),
                "n_points": int(stats["n_points"]),
                "n_significant_points": int(stats["n_sig"]),
                "n_recurrent_gene_hits": int(stats["n_recurrent_hits"]),
            }
        )
    volcano_manifest = pd.DataFrame(volcano_rows)

    ora_terms_selected, ora_summary_selected = _subset_ora_for_targets(
        ora_tables=args.ora_tables,
        ora_streams=ora_streams,
        targets=targets_present,
    )

    chunk_summary_df = (
        pd.read_csv(args.chunk_summary, sep="\t", compression="infer")
        if args.chunk_summary.exists()
        else pd.DataFrame()
    )

    selected_out = args.outdir / "selected_perturbations.tsv"
    chunk_summary_out = args.outdir / "chunk_summary.tsv"
    expr_out = args.outdir / "pseudobulk_expression_matrix.tsv.gz"
    deg_out = args.outdir / "deg_profile_matrix.tsv.gz"
    expr_clusters_out = args.outdir / "perturbation_clusters_expression.tsv"
    deg_clusters_out = args.outdir / "perturbation_clusters.tsv"
    expr_heatmap_out = args.outdir / "heatmap_expression.png"
    deg_heatmap_out = args.outdir / "heatmap_deg_profile.png"
    volcano_manifest_out = args.outdir / "volcano_plots.tsv"
    ora_terms_targets_out = args.outdir / "ora_terms_selected_targets.tsv.gz"
    ora_summary_targets_out = args.outdir / "ora_summary_selected_targets.tsv.gz"
    meta_out = args.outdir / "clustering_meta.json"
    done_out = args.outdir / "done.txt"
    methods_summary_out = args.outdir / "method_clustering_summary.tsv"

    selected_df.to_csv(selected_out, sep="\t", index=False)
    chunk_summary_df.to_csv(chunk_summary_out, sep="\t", index=False)
    expr_mat.reset_index().to_csv(expr_out, sep="\t", index=False, compression="gzip")
    deg_mat.reset_index().to_csv(deg_out, sep="\t", index=False, compression="gzip")
    expr_clusters.to_csv(expr_clusters_out, sep="\t", index=False)
    deg_clusters.to_csv(deg_clusters_out, sep="\t", index=False)
    volcano_manifest.to_csv(volcano_manifest_out, sep="\t", index=False)
    ora_terms_selected.to_csv(ora_terms_targets_out, sep="\t", index=False, compression="gzip")
    ora_summary_selected.to_csv(ora_summary_targets_out, sep="\t", index=False, compression="gzip")
    _write_pairwise_heatmap(
        matrix=expr_pair,
        order=expr_ordered,
        groups=row_groups,
        out_png=expr_heatmap_out,
        title=f"{args.dataset}: Pseudobulk expression similarity (rows/cols marked)",
    )
    _write_pairwise_heatmap(
        matrix=deg_pair,
        order=deg_ordered,
        groups=row_groups,
        out_png=deg_heatmap_out,
        title=f"{args.dataset}: DEG-profile similarity (rows/cols marked)",
    )

    methods_root = args.outdir / "methods"
    methods_root.mkdir(parents=True, exist_ok=True)
    method_summary_rows: list[dict[str, object]] = []

    for expr_method_name, expr_assignment_source in [("mixscape", mixscape_map), ("mixscale", mixscale_map)]:
        method_label = f"expression_{expr_method_name}"
        method_dir = methods_root / method_label
        method_dir.mkdir(parents=True, exist_ok=True)
        rng_method = np.random.default_rng(_seed_for_label(int(args.random_seed), method_label))
        targets_present_m, missing_targets_m, sampled_random_m = _sample_target_and_random(
            requested_targets=requested_targets,
            candidate_perts=set(expr_assignment_source.keys()),
            random_match_target_count=bool(args.random_match_target_count),
            n_random_controls=int(args.n_random_controls),
            rng=rng_method,
        )
        selected_ko_m = list(targets_present_m) + list(sampled_random_m)
        assignment_map_m = {p: set(expr_assignment_source.get(p, set())) for p in selected_ko_m if expr_assignment_source.get(p, set())}
        selected_ko_m = [p for p in selected_ko_m if p in assignment_map_m]

        expr_mat_m, expr_used_counts_m = _build_pseudobulk_expression(
            h5ad=args.h5ad,
            assignment_map=assignment_map_m,
            control_label=str(args.control_label),
            control_barcodes=control_barcodes,
            normalize_target_sum=float(args.expression_normalize_target_sum),
            log1p=bool(args.expression_log1p),
        )
        selected_df_m = _build_selected_df(
            selected_ko_perts=selected_ko_m,
            targets_present=targets_present_m,
            sampled_random=sampled_random_m,
            control_label=str(args.control_label),
            present_in_expression=set(expr_mat_m.index.astype(str).tolist()),
            present_in_de=set(),
            used_counts=expr_used_counts_m,
        )
        clusters_m, ordered_m = _cluster_matrix(
            matrix=expr_mat_m,
            n_clusters=int(args.n_clusters),
            linkage_method=str(args.linkage),
            distance_metric=str(args.distance),
        )
        clusters_m = clusters_m.merge(selected_df_m[["perturbation", "group", "is_target", "is_control"]], on="perturbation", how="left")
        pair_m = _pairwise_similarity(expr_mat_m)
        heatmap_mat_m = _select_heatmap_columns(
            matrix=expr_mat_m,
            max_genes=int(args.heatmap_max_genes),
            gene_order=str(args.heatmap_gene_order),
        )
        row_groups_m = dict(zip(selected_df_m["perturbation"].astype(str), selected_df_m["group"].astype(str)))

        selected_out_m = method_dir / "selected_perturbations.tsv"
        matrix_out_m = method_dir / "pseudobulk_expression_matrix.tsv.gz"
        clusters_out_m = method_dir / "perturbation_clusters.tsv"
        heatmap_out_m = method_dir / "heatmap_expression.png"
        meta_out_m = method_dir / "clustering_meta.json"

        selected_df_m.to_csv(selected_out_m, sep="\t", index=False)
        expr_mat_m.reset_index().to_csv(matrix_out_m, sep="\t", index=False, compression="gzip")
        clusters_m.to_csv(clusters_out_m, sep="\t", index=False)
        _write_pairwise_heatmap(
            matrix=pair_m,
            order=ordered_m,
            groups=row_groups_m,
            out_png=heatmap_out_m,
            title=f"{args.dataset}: Pseudobulk expression similarity ({expr_method_name})",
        )
        method_meta_m = {
            "dataset": args.dataset,
            "method": method_label,
            "mode": "expression",
            "expression_source": expr_method_name,
            "target_list": str(args.target_list),
            "n_targets_requested": int(len(requested_targets)),
            "n_targets_present": int(len(targets_present_m)),
            "n_random_controls_sampled": int(len(sampled_random_m)),
            "n_selected_ko_perturbations": int(len(selected_ko_m)),
            "n_selected_rows_including_control": int(selected_df_m.shape[0]),
            "missing_targets": missing_targets_m,
            "n_control_cells_used": int(expr_used_counts_m.get(str(args.control_label), 0)),
            "ko_cells_used_per_perturbation": {k: int(v) for k, v in expr_used_counts_m.items() if k != str(args.control_label)},
            "n_genes_expression_profile": int(expr_mat_m.shape[1]),
            "n_genes_expression_heatmap": int(heatmap_mat_m.shape[1]),
            "ordered_perturbations_expression": ordered_m,
        }
        meta_out_m.write_text(json.dumps(method_meta_m, indent=2))
        method_summary_rows.append(
            {
                "method": method_label,
                "mode": "expression",
                "n_targets_requested": int(len(requested_targets)),
                "n_targets_present": int(len(targets_present_m)),
                "n_random_controls_sampled": int(len(sampled_random_m)),
                "n_selected_ko_perturbations": int(len(selected_ko_m)),
                "n_genes_profile": int(expr_mat_m.shape[1]),
                "selected_path": str(selected_out_m),
                "matrix_path": str(matrix_out_m),
                "clusters_path": str(clusters_out_m),
                "heatmap_path": str(heatmap_out_m),
                "meta_path": str(meta_out_m),
            }
        )

    for deg_method_name in ["global_de", "mixscape", "mixscale"]:
        method_label = f"deg_{deg_method_name}"
        method_dir = methods_root / method_label
        method_dir.mkdir(parents=True, exist_ok=True)
        rng_method = np.random.default_rng(_seed_for_label(int(args.random_seed), method_label))
        de_agg_m = de_agg_by_stream.get(deg_method_name, _aggregate_standardized_deg(pd.DataFrame(), denom=1))
        candidate_perts_m = set(de_agg_m["perturbation"].astype(str).tolist()) if not de_agg_m.empty else set()
        targets_present_m, missing_targets_m, sampled_random_m = _sample_target_and_random(
            requested_targets=requested_targets,
            candidate_perts=candidate_perts_m,
            random_match_target_count=bool(args.random_match_target_count),
            n_random_controls=int(args.n_random_controls),
            rng=rng_method,
        )
        selected_ko_m = list(targets_present_m) + list(sampled_random_m)
        deg_mat_m = _build_deg_matrix(
            de_agg=de_agg_m,
            selected_perts=selected_ko_m,
            control_label=str(args.control_label),
            top_genes_per_perturbation=int(args.top_genes_per_perturbation),
            min_feature_genes=int(args.min_feature_genes),
        )
        selected_df_m = _build_selected_df(
            selected_ko_perts=selected_ko_m,
            targets_present=targets_present_m,
            sampled_random=sampled_random_m,
            control_label=str(args.control_label),
            present_in_expression=set(),
            present_in_de=set(deg_mat_m.index.astype(str).tolist()),
            used_counts={},
        )
        clusters_m, ordered_m = _cluster_matrix(
            matrix=deg_mat_m,
            n_clusters=int(args.n_clusters),
            linkage_method=str(args.linkage),
            distance_metric=str(args.distance),
        )
        clusters_m = clusters_m.merge(selected_df_m[["perturbation", "group", "is_target", "is_control"]], on="perturbation", how="left")
        pair_m = _pairwise_similarity(deg_mat_m)
        heatmap_mat_m = _select_heatmap_columns(
            matrix=deg_mat_m,
            max_genes=int(args.heatmap_max_genes),
            gene_order=str(args.heatmap_gene_order),
        )
        row_groups_m = dict(zip(selected_df_m["perturbation"].astype(str), selected_df_m["group"].astype(str)))

        selected_out_m = method_dir / "selected_perturbations.tsv"
        matrix_out_m = method_dir / "deg_profile_matrix.tsv.gz"
        clusters_out_m = method_dir / "perturbation_clusters.tsv"
        heatmap_out_m = method_dir / "heatmap_deg_profile.png"
        meta_out_m = method_dir / "clustering_meta.json"

        selected_df_m.to_csv(selected_out_m, sep="\t", index=False)
        deg_mat_m.reset_index().to_csv(matrix_out_m, sep="\t", index=False, compression="gzip")
        clusters_m.to_csv(clusters_out_m, sep="\t", index=False)
        _write_pairwise_heatmap(
            matrix=pair_m,
            order=ordered_m,
            groups=row_groups_m,
            out_png=heatmap_out_m,
            title=f"{args.dataset}: DEG-profile similarity ({deg_method_name})",
        )
        method_meta_m = {
            "dataset": args.dataset,
            "method": method_label,
            "mode": "deg",
            "deg_stream": deg_method_name,
            "stream_available": bool(deg_method_name in stream_names),
            "target_list": str(args.target_list),
            "n_targets_requested": int(len(requested_targets)),
            "n_targets_present": int(len(targets_present_m)),
            "n_random_controls_sampled": int(len(sampled_random_m)),
            "n_selected_ko_perturbations": int(len(selected_ko_m)),
            "n_selected_rows_including_control": int(selected_df_m.shape[0]),
            "missing_targets": missing_targets_m,
            "n_genes_deg_profile": int(deg_mat_m.shape[1]),
            "n_genes_deg_heatmap": int(heatmap_mat_m.shape[1]),
            "ordered_perturbations_deg": ordered_m,
        }
        meta_out_m.write_text(json.dumps(method_meta_m, indent=2))
        method_summary_rows.append(
            {
                "method": method_label,
                "mode": "deg",
                "n_targets_requested": int(len(requested_targets)),
                "n_targets_present": int(len(targets_present_m)),
                "n_random_controls_sampled": int(len(sampled_random_m)),
                "n_selected_ko_perturbations": int(len(selected_ko_m)),
                "n_genes_profile": int(deg_mat_m.shape[1]),
                "selected_path": str(selected_out_m),
                "matrix_path": str(matrix_out_m),
                "clusters_path": str(clusters_out_m),
                "heatmap_path": str(heatmap_out_m),
                "meta_path": str(meta_out_m),
            }
        )

    method_summary_df = pd.DataFrame(method_summary_rows)
    method_summary_df.to_csv(methods_summary_out, sep="\t", index=False)

    meta = {
        "dataset": args.dataset,
        "deg_streams": stream_names,
        "de_paths": [str(p) for p in args.de_tables],
        "de_stream_meta": de_stream_meta,
        "n_streams_total": int(n_streams_total),
        "ora_streams_linked": ora_streams,
        "ora_paths": [str(p) for p in args.ora_tables],
        "mixscape_selected_path": str(args.mixscape_selected),
        "mixscale_cells_path": str(args.mixscale_cells),
        "assignment_logic": str(args.assignment_logic),
        "mixscale_score_threshold": float(args.mixscale_score_threshold),
        "h5ad": str(args.h5ad),
        "pert_col": args.pert_col,
        "control_label": str(args.control_label),
        "target_list": str(args.target_list),
        "n_targets_requested": int(len(requested_targets)),
        "n_targets_present": int(len(targets_present)),
        "n_random_controls_sampled": int(len(sampled_random)),
        "n_selected_ko_perturbations": int(len(selected_ko_perts)),
        "n_selected_rows_including_control": int(selected_df.shape[0]),
        "n_genes_expression_profile": int(expr_mat.shape[1]),
        "n_genes_deg_profile": int(deg_mat.shape[1]),
        "n_genes_expression_heatmap": int(expr_heatmap_mat.shape[1]),
        "n_genes_deg_heatmap": int(deg_heatmap_mat.shape[1]),
        "ordered_perturbations_expression": expr_ordered,
        "ordered_perturbations_deg": deg_ordered,
        "missing_targets": missing_targets,
        "n_control_cells_used": int(expr_used_counts.get(str(args.control_label), 0)),
        "ko_cells_used_per_perturbation": {k: int(v) for k, v in expr_used_counts.items() if k != str(args.control_label)},
        "volcano_targets_plotted": volc_targets,
        "recurrent_gene_fraction": float(args.recurrent_gene_fraction),
        "recurrent_gene_count": int(len(recurrent_genes)),
        "ora_terms_rows_selected_targets": int(ora_terms_selected.shape[0]),
        "ora_summary_rows_selected_targets": int(ora_summary_selected.shape[0]),
        "filters": {
            "random_controls_match_target_count": bool(args.random_match_target_count),
            "n_random_controls": int(args.n_random_controls),
            "max_pval_adj": float(args.max_pval_adj),
            "min_abs_effect": float(args.min_abs_effect),
            "rank_significant_max": int(args.rank_significant_max),
            "top_genes_per_perturbation": int(args.top_genes_per_perturbation),
            "min_feature_genes": int(args.min_feature_genes),
            "max_control_cells": int(args.max_control_cells),
            "expression_normalize_target_sum": float(args.expression_normalize_target_sum),
            "expression_log1p": bool(args.expression_log1p),
            "heatmap_max_genes": int(args.heatmap_max_genes),
            "heatmap_gene_order": str(args.heatmap_gene_order),
            "volcano_max_targets": int(args.volcano_max_targets),
            "volcano_stream_priority": volc_stream_priority,
            "n_clusters": int(args.n_clusters),
            "linkage": str(args.linkage),
            "distance": str(args.distance),
            "random_seed": int(args.random_seed),
        },
        "n_de_rows_raw_total": int(de_raw_all.shape[0]),
        "n_de_rows_aggregated": int(de_agg.shape[0]),
        "n_chunks_from_summary": int(chunk_summary_df.shape[0]),
        "method_outputs_root": str(methods_root),
        "method_clustering_summary": str(methods_summary_out),
        "method_clustering_rows": int(method_summary_df.shape[0]),
    }
    meta_out.write_text(json.dumps(meta, indent=2))
    done_out.write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", type=str, default="dataset")
    ap.add_argument("--de-tables", nargs="+", type=Path, required=True)
    ap.add_argument("--deg-streams", type=str, default="global_de,mixscape,mixscale")
    ap.add_argument("--ora-tables", nargs="*", type=Path, default=[])
    ap.add_argument("--ora-streams", type=str, default="")
    ap.add_argument("--mixscape-selected", type=Path, required=True)
    ap.add_argument("--mixscale-cells", type=Path, required=True)
    ap.add_argument("--chunk-summary", type=Path, required=True)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--target-list", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--random-match-target-count", action="store_true")
    ap.add_argument("--n-random-controls", type=int, default=0)
    ap.add_argument("--min-abs-effect", type=float, default=0.25)
    ap.add_argument("--max-pval-adj", type=float, default=0.05)
    ap.add_argument("--rank-significant-max", type=int, default=100)
    ap.add_argument("--top-genes-per-perturbation", type=int, default=0)
    ap.add_argument("--min-feature-genes", type=int, default=20)
    ap.add_argument("--assignment-logic", type=str, default="intersection")
    ap.add_argument("--mixscale-score-threshold", type=float, default=0.0)
    ap.add_argument("--max-control-cells", type=int, default=50000)
    ap.add_argument("--expression-normalize-target-sum", type=float, default=10000)
    ap.add_argument("--expression-log1p", dest="expression_log1p", action="store_true")
    ap.add_argument("--no-expression-log1p", dest="expression_log1p", action="store_false")
    ap.set_defaults(expression_log1p=True)
    ap.add_argument("--heatmap-max-genes", type=int, default=200)
    ap.add_argument("--heatmap-gene-order", type=str, default="variance")
    ap.add_argument("--volcano-max-targets", type=int, default=20)
    ap.add_argument("--recurrent-gene-fraction", type=float, default=0.25)
    ap.add_argument("--volcano-stream-priority", type=str, default="global_de,mixscape,ps")
    ap.add_argument("--n-clusters", type=int, default=20)
    ap.add_argument("--linkage", type=str, default="average")
    ap.add_argument("--distance", type=str, default="euclidean")
    ap.add_argument("--random-seed", type=int, default=0)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        dataset=str(snk.params.dataset),
        de_tables=[Path(str(p)) for p in list(snk.input.de_tables)],
        deg_streams=str(snk.params.deg_streams),
        ora_tables=[Path(str(p)) for p in list(snk.input.ora_tables)],
        ora_streams=str(snk.params.ora_streams),
        mixscape_selected=Path(str(snk.input.mixscape_selected)),
        mixscale_cells=Path(str(snk.input.mixscale_cells)),
        chunk_summary=Path(str(snk.input.chunk_summary)),
        h5ad=Path(str(snk.input.h5ad)),
        target_list=Path(str(snk.input.target_list)),
        outdir=Path(str(snk.params.outdir)),
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        random_match_target_count=bool(snk.params.random_match_target_count),
        n_random_controls=int(snk.params.n_random),
        min_abs_effect=float(snk.params.min_abs_effect),
        max_pval_adj=float(snk.params.max_pval_adj),
        rank_significant_max=int(snk.params.rank_significant_max),
        top_genes_per_perturbation=int(snk.params.top_genes_per_pert),
        min_feature_genes=int(snk.params.min_feature_genes),
        assignment_logic=str(snk.params.assignment_logic),
        mixscale_score_threshold=float(snk.params.mixscale_score_threshold),
        max_control_cells=int(snk.params.max_control_cells),
        expression_normalize_target_sum=float(snk.params.expression_normalize_target_sum),
        expression_log1p=bool(snk.params.expression_log1p),
        heatmap_max_genes=int(snk.params.heatmap_max_genes),
        heatmap_gene_order=str(snk.params.heatmap_gene_order),
        volcano_max_targets=int(snk.params.volcano_max_targets),
        recurrent_gene_fraction=float(snk.params.recurrent_gene_fraction),
        volcano_stream_priority=str(snk.params.volcano_stream_priority),
        n_clusters=int(snk.params.n_clusters),
        linkage=str(snk.params.linkage),
        distance=str(snk.params.distance),
        random_seed=int(snk.params.seed),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
