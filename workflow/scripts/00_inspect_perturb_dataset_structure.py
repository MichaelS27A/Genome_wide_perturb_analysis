#!/usr/bin/env python3
"""Inspect perturbation dataset structure for Snakemake Mixscape pipeline readiness.

This script is intended for running in external folders.
It checks .h5ad structure and optional guide-calls CSV compatibility.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


def summarize_h5ad(h5ad_path: Path, pert_col: str, control_label: str) -> dict[str, Any]:
    # Try scanpy-backed first. Some very large h5ad layouts can still trigger
    # heavy layer materialization in some anndata versions, so we fall back to
    # direct HDF5 inspection if needed.
    try:
        import scanpy as sc

        adata = sc.read_h5ad(h5ad_path, backed="r")
        try:
            obs = adata.obs
            var = adata.var

            info: dict[str, Any] = {
                "path": str(h5ad_path),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "obs_columns": [str(c) for c in obs.columns.tolist()],
                "var_columns": [str(c) for c in var.columns.tolist()],
                "layers": sorted(list(adata.layers.keys())),
                "obsm": sorted(list(adata.obsm.keys())),
                "uns_keys": sorted(list(adata.uns.keys())),
                "pert_col": pert_col,
                "pert_col_present": bool(pert_col in obs.columns),
                "read_mode": "scanpy_backed",
            }

            if pert_col in obs.columns:
                ser = obs[pert_col].astype("string")
                vc = ser.value_counts(dropna=False).head(20)
                info["top_perturbations"] = [{"label": str(k), "count": int(v)} for k, v in vc.items()]
                info["n_unique_perturbations"] = int(ser.nunique(dropna=True))
                info["control_label"] = control_label
                info["control_label_count"] = int((ser == control_label).sum())
            else:
                info["top_perturbations"] = []
                info["n_unique_perturbations"] = 0
                info["control_label"] = control_label
                info["control_label_count"] = 0

            return info
        finally:
            adata.file.close()
    except Exception as e:
        return summarize_h5ad_via_h5py(h5ad_path, pert_col, control_label, fallback_error=str(e))


def summarize_h5ad_via_h5py(
    h5ad_path: Path, pert_col: str, control_label: str, fallback_error: str
) -> dict[str, Any]:
    with h5py.File(h5ad_path, "r") as f:
        obs_keys = sorted(list(f["obs"].keys())) if "obs" in f else []
        var_keys = sorted(list(f["var"].keys())) if "var" in f else []
        layers = sorted(list(f["layers"].keys())) if "layers" in f else []
        obsm = sorted(list(f["obsm"].keys())) if "obsm" in f else []
        uns_keys = sorted(list(f["uns"].keys())) if "uns" in f else []

        n_obs = len(f["obs"]["_index"]) if "obs" in f and "_index" in f["obs"] else 0
        n_vars = len(f["var"]["_index"]) if "var" in f and "_index" in f["var"] else 0

        info: dict[str, Any] = {
            "path": str(h5ad_path),
            "n_obs": int(n_obs),
            "n_vars": int(n_vars),
            "obs_columns": [str(c) for c in obs_keys if c != "_index"],
            "var_columns": [str(c) for c in var_keys if c != "_index"],
            "layers": layers,
            "obsm": obsm,
            "uns_keys": uns_keys,
            "pert_col": pert_col,
            "pert_col_present": bool(pert_col in obs_keys),
            "read_mode": "h5py_fallback",
            "fallback_error": fallback_error,
            "top_perturbations": [],
            "n_unique_perturbations": 0,
            "control_label": control_label,
            "control_label_count": 0,
        }

        if pert_col in obs_keys:
            obj = f["obs"][pert_col]
            if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
                categories = obj["categories"][:]
                categories = [
                    c.decode() if isinstance(c, (bytes, np.bytes_)) else str(c) for c in categories
                ]
                codes = obj["codes"][:]
                valid = codes[codes >= 0]
                vals, counts = np.unique(valid, return_counts=True)
                order = np.argsort(counts)[::-1]

                top = []
                for idx in order[:20]:
                    cat_idx = int(vals[idx])
                    top.append({"label": categories[cat_idx], "count": int(counts[idx])})
                info["top_perturbations"] = top
                info["n_unique_perturbations"] = int(len(categories))

                if control_label in categories:
                    cidx = categories.index(control_label)
                    ccount = int(counts[vals.tolist().index(cidx)]) if cidx in vals else 0
                    info["control_label_count"] = ccount
            elif isinstance(obj, h5py.Dataset):
                arr = obj[:]
                ser = pd.Series(arr).astype("string")
                vc = ser.value_counts(dropna=False).head(20)
                info["top_perturbations"] = [{"label": str(k), "count": int(v)} for k, v in vc.items()]
                info["n_unique_perturbations"] = int(ser.nunique(dropna=True))
                info["control_label_count"] = int((ser == control_label).sum())

        return info


def summarize_guide_csv(
    csv_path: Path,
    obs_names: pd.Index,
    control_label: str,
    max_rows: int,
    chunk_size: int,
) -> dict[str, Any]:
    required = {"cell_barcode", "gene_target"}
    out: dict[str, Any] = {
        "path": str(csv_path),
        "required_columns": sorted(required),
        "columns_present": [],
        "missing_columns": [],
        "rows_scanned": 0,
        "unique_barcodes_scanned": 0,
        "unique_gene_targets_scanned": 0,
        "control_label": control_label,
        "control_rows_scanned": 0,
        "obs_overlap_barcodes_scanned": 0,
    }

    barcode_set: set[str] = set()
    gene_counts: dict[str, int] = {}

    reader = pd.read_csv(csv_path, compression="infer", chunksize=chunk_size, low_memory=True)
    for chunk in reader:
        if out["rows_scanned"] >= max_rows:
            break

        if not out["columns_present"]:
            out["columns_present"] = [str(c) for c in chunk.columns.tolist()]
            missing = sorted(list(required - set(chunk.columns)))
            out["missing_columns"] = missing
            if missing:
                return out

        need = min(len(chunk), max_rows - out["rows_scanned"])
        ch = chunk.iloc[:need][["cell_barcode", "gene_target"]].copy()

        ch["cell_barcode"] = ch["cell_barcode"].astype("string")
        ch["gene_target"] = ch["gene_target"].astype("string")

        out["rows_scanned"] += int(ch.shape[0])
        out["control_rows_scanned"] += int((ch["gene_target"] == control_label).sum())

        bcs = ch["cell_barcode"].dropna().astype(str).tolist()
        barcode_set.update(bcs)

        vc = ch["gene_target"].fillna("<NA>").astype(str).value_counts(dropna=False)
        for k, v in vc.items():
            gene_counts[str(k)] = gene_counts.get(str(k), 0) + int(v)

    out["unique_barcodes_scanned"] = len(barcode_set)
    out["unique_gene_targets_scanned"] = len(gene_counts)
    out["obs_overlap_barcodes_scanned"] = int(pd.Index(list(barcode_set)).isin(obs_names).sum())

    top_targets = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:20]
    out["top_gene_targets_scanned"] = [{"label": str(k), "count": int(v)} for k, v in top_targets]
    return out


def build_readiness_report(h5ad: dict[str, Any], guide: dict[str, Any] | None) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    checks.append(
        {
            "name": "h5ad_pert_col_present",
            "ok": bool(h5ad["pert_col_present"]),
            "detail": f"pert_col={h5ad['pert_col']}",
        }
    )
    checks.append(
        {
            "name": "h5ad_has_controls",
            "ok": bool(h5ad["control_label_count"] > 0),
            "detail": f"control_label={h5ad['control_label']} count={h5ad['control_label_count']}",
        }
    )

    if guide is not None:
        checks.append(
            {
                "name": "guide_csv_required_columns",
                "ok": len(guide.get("missing_columns", [])) == 0,
                "detail": f"missing={guide.get('missing_columns', [])}",
            }
        )
        checks.append(
            {
                "name": "guide_csv_obs_barcode_overlap",
                "ok": bool(guide.get("obs_overlap_barcodes_scanned", 0) > 0),
                "detail": f"overlap_scanned={guide.get('obs_overlap_barcodes_scanned', 0)}",
            }
        )

    ok_all = all(bool(c["ok"]) for c in checks)
    return {"ready": ok_all, "checks": checks}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--guide-calls-csv", type=Path, default=None)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--max-csv-rows", type=int, default=1_000_000)
    ap.add_argument("--csv-chunk-size", type=int, default=200_000)
    ap.add_argument("--output-json", type=Path, default=None)
    args = ap.parse_args()

    if not args.h5ad.exists():
        raise FileNotFoundError(f"h5ad not found: {args.h5ad}")
    if args.guide_calls_csv and not args.guide_calls_csv.exists():
        raise FileNotFoundError(f"guide-calls CSV not found: {args.guide_calls_csv}")

    h5ad_summary = summarize_h5ad(args.h5ad, args.pert_col, args.control_label)
    obs_names = pd.Index([])

    # We only need obs names for overlap checks.
    try:
        import scanpy as sc

        adata = sc.read_h5ad(args.h5ad, backed="r")
        try:
            obs_names = pd.Index(adata.obs_names.astype(str))
        finally:
            adata.file.close()
    except Exception:
        with h5py.File(args.h5ad, "r") as f:
            if "obs" in f and "_index" in f["obs"]:
                idx = f["obs"]["_index"][:]
                obs_names = pd.Index(
                    [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in idx]
                )

    guide_summary = None
    if args.guide_calls_csv is not None:
        guide_summary = summarize_guide_csv(
            csv_path=args.guide_calls_csv,
            obs_names=obs_names,
            control_label=args.control_label,
            max_rows=args.max_csv_rows,
            chunk_size=args.csv_chunk_size,
        )

    readiness = build_readiness_report(h5ad_summary, guide_summary)
    result = {
        "h5ad": h5ad_summary,
        "guide_calls_csv": guide_summary,
        "pipeline_readiness": readiness,
    }

    if args.output_json is not None:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(json.dumps(result, indent=2))

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
