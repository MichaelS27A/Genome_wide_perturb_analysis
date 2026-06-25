#!/usr/bin/env python3
"""Compare Mixscape, Mixscale and PS outputs for one dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer")


def _z(x: pd.Series) -> pd.Series:
    if x.empty:
        return x
    sd = float(x.std(ddof=0))
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(x)), index=x.index)
    return (x - float(x.mean())) / sd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--mixscape-stats", type=Path, required=True)
    ap.add_argument("--mixscale-cell-scores", type=Path, required=True)
    ap.add_argument("--mixscale-de", type=Path, required=True)
    ap.add_argument("--ps-cell-scores", type=Path, required=True)
    ap.add_argument("--ps-summary", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    ms = _read_tsv(args.mixscape_stats)
    msc_cell = _read_tsv(args.mixscale_cell_scores)
    msc_de = _read_tsv(args.mixscale_de)
    ps_cell = _read_tsv(args.ps_cell_scores)
    ps_sum = _read_tsv(args.ps_summary)

    msc_effect_col = "effect_l2_vs_control" if "effect_l2_vs_control" in ms.columns else None
    if msc_effect_col is None:
        raise RuntimeError("mixscape stats missing expected column 'effect_l2_vs_control'")

    mixscape_tbl = ms[["perturbation", msc_effect_col]].rename(columns={msc_effect_col: "mixscape_effect"})

    mixscale_cell_agg = (
        msc_cell.groupby("perturbation", as_index=False)["mixscale_score"]
        .mean()
        .rename(columns={"mixscale_score": "mixscale_score_mean"})
    )

    if "p_weight" in msc_de.columns:
        mixscale_de_agg = (
            msc_de.groupby("perturbation", as_index=False)["p_weight"]
            .min()
            .rename(columns={"p_weight": "mixscale_best_p"})
        )
    else:
        mixscale_de_agg = pd.DataFrame({"perturbation": mixscale_cell_agg["perturbation"], "mixscale_best_p": np.nan})

    ps_cell_agg = (
        ps_cell.groupby("perturbation", as_index=False)["ps_score"]
        .mean()
        .rename(columns={"ps_score": "ps_score_mean"})
    )

    long = (
        mixscape_tbl.merge(mixscale_cell_agg, on="perturbation", how="outer")
        .merge(mixscale_de_agg, on="perturbation", how="outer")
        .merge(ps_cell_agg, on="perturbation", how="outer")
        .merge(ps_sum[["perturbation", "mean_ps_score"]], on="perturbation", how="left")
    )

    for col in ["mixscape_effect", "mixscale_score_mean", "ps_score_mean", "mean_ps_score"]:
        if col in long.columns:
            long[f"{col}_z"] = _z(pd.to_numeric(long[col], errors="coerce"))

    corr_rows = []
    metric_cols = [c for c in ["mixscape_effect", "mixscale_score_mean", "ps_score_mean", "mean_ps_score"] if c in long.columns]
    for i, a in enumerate(metric_cols):
        for b in metric_cols[i + 1 :]:
            ab = long[[a, b]].dropna()
            corr = float(ab[a].corr(ab[b])) if len(ab) >= 3 else np.nan
            corr_rows.append({"metric_a": a, "metric_b": b, "pearson_r": corr, "n": int(len(ab))})

    summary = pd.DataFrame(corr_rows)

    long.to_csv(args.outdir / "comparison_long.tsv.gz", sep="\t", index=False, compression="gzip")
    summary.to_csv(args.outdir / "comparison_summary.tsv", sep="\t", index=False)

    meta = {
        "dataset": args.dataset,
        "n_perturbations": int(long["perturbation"].nunique(dropna=True)),
        "inputs": {
            "mixscape_stats": str(args.mixscape_stats),
            "mixscale_cell_scores": str(args.mixscale_cell_scores),
            "mixscale_de": str(args.mixscale_de),
            "ps_cell_scores": str(args.ps_cell_scores),
            "ps_summary": str(args.ps_summary),
        },
    }
    (args.outdir / "method_meta.json").write_text(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
