#!/usr/bin/env python3
"""Batch inspect .h5ad dataset structure and suggest pert_col/control_label.

This script is intentionally lightweight and h5py-first so it can inspect very
large .h5ad files without loading expression matrices into memory.
"""


import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd

try:
    import yaml
except Exception:  # pragma: no cover
    yaml = None


DEFAULT_DATASETS: dict[str, str] = {
    "HCT116": "data/HCT116_filtered_dual_guide_cells.h5ad",
    "HEK293T": "data/HEK293T_filtered_dual_guide_cells.h5ad",
    "STEMCELL": "data/STEMCELL.h5ad",
    "K562_GWPS_RAW_FULL": "data/K562_GWPS_RAW_FULL.h5ad",
}

CANDIDATE_PERT_COLS = [
    "gene_target",
    "guide_target",
    "target_gene",
    "perturbation",
    "perturbation_label",
    "gene",
    "target",
    "sgRNA",
    "sgrna",
    "guide",
]

CANDIDATE_CONTROL_LABELS = [
    "Non-Targeting",
    "non-targeting",
    "NTC",
    "NT",
    "negative",
    "neg",
    "control",
    "CTRL",
]


@dataclass
class DatasetSpec:
    name: str
    h5ad: Path


def _decode_value(x: Any) -> str:
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode()
    return str(x)


def _top_counts_from_codes(categories: list[str], codes: np.ndarray, top_n: int = 20) -> list[dict[str, Any]]:
    valid = codes[codes >= 0]
    if valid.size == 0:
        return []
    vals, counts = np.unique(valid, return_counts=True)
    order = np.argsort(counts)[::-1][:top_n]
    out: list[dict[str, Any]] = []
    for idx in order:
        cat_idx = int(vals[idx])
        out.append({"label": categories[cat_idx], "count": int(counts[idx])})
    return out


def inspect_obs_column(
    f: h5py.File,
    col: str,
    control_labels: list[str],
    max_plain_rows: int = 2_000_000,
) -> dict[str, Any]:
    if "obs" not in f or col not in f["obs"]:
        return {"present": False, "name": col}

    obj = f["obs"][col]
    result: dict[str, Any] = {
        "present": True,
        "name": col,
        "encoding": "unknown",
        "n_unique": 0,
        "top_values": [],
        "control_counts": {lbl: 0 for lbl in control_labels},
    }

    # AnnData categorical format: obs/<col> as group with categories + codes.
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        result["encoding"] = "categorical_codes"
        categories = [_decode_value(x) for x in obj["categories"][:]]
        codes = obj["codes"][:]
        result["n_unique"] = int(len(categories))
        result["top_values"] = _top_counts_from_codes(categories, codes)

        if categories:
            for lbl in control_labels:
                if lbl in categories:
                    idx = categories.index(lbl)
                    result["control_counts"][lbl] = int(np.sum(codes == idx))
        return result

    # Plain dataset fallback; sample if extremely large.
    if isinstance(obj, h5py.Dataset):
        result["encoding"] = "dataset"
        n = int(obj.shape[0]) if obj.shape else 0
        take = min(n, max_plain_rows)
        arr = obj[:take]

        # Legacy AnnData (0.1.0 dataframe encoding) stores category labels in
        # obs/__categories/<col> and integer codes in obs/<col>.
        legacy_categories = None
        if "obs" in f and "__categories" in f["obs"]:
            cats_group = f["obs"]["__categories"]
            if isinstance(cats_group, h5py.Group) and col in cats_group:
                try:
                    legacy_categories = [_decode_value(x) for x in cats_group[col][:]]
                except Exception:
                    legacy_categories = None

        if legacy_categories is not None and np.issubdtype(arr.dtype, np.integer):
            ser = pd.Series(
                [
                    legacy_categories[int(i)] if 0 <= int(i) < len(legacy_categories) else "<NA>"
                    for i in arr
                ],
                dtype="string",
            )
            result["encoding"] = "legacy_categorical_codes"
        else:
            ser = pd.Series([_decode_value(x) for x in arr], dtype="string")

        vc = ser.value_counts(dropna=False).head(20)
        result["n_unique"] = int(ser.nunique(dropna=True))
        result["top_values"] = [{"label": str(k), "count": int(v)} for k, v in vc.items()]
        for lbl in control_labels:
            result["control_counts"][lbl] = int((ser == lbl).sum())
        result["sampled_rows"] = int(take)
        result["total_rows"] = int(n)
        return result

    return result


def _sample_numeric_from_dataset(ds: h5py.Dataset, sample_n: int = 200_000) -> np.ndarray:
    if not ds.shape:
        return np.array([], dtype=float)

    # 1D vectors (e.g., sparse matrix data arrays)
    if len(ds.shape) == 1:
        n = int(ds.shape[0])
        if n <= 0:
            return np.array([], dtype=float)
        take = min(n, sample_n)
        arr = ds[:take]
    else:
        # 2D+ arrays: sample a bounded top-left block to avoid huge allocations.
        n0 = int(ds.shape[0])
        n1 = int(ds.shape[1])
        if n0 <= 0 or n1 <= 0:
            return np.array([], dtype=float)
        max_rows = min(n0, 1000)
        max_cols = min(n1, max(1, sample_n // max_rows))
        arr = ds[:max_rows, :max_cols]

    try:
        flat = np.asarray(arr, dtype=float).reshape(-1)
        if flat.size > sample_n:
            return flat[:sample_n]
        return flat
    except Exception:
        return np.array([], dtype=float)


def _integer_like_fraction(arr: np.ndarray) -> float | None:
    if arr.size == 0:
        return None
    if not np.issubdtype(arr.dtype, np.number):
        return None
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return None
    return float(np.mean(np.isclose(finite, np.round(finite), atol=1e-8)))


def _inspect_matrix_obj(obj: Any) -> dict[str, Any]:
    out: dict[str, Any] = {"kind": type(obj).__name__}
    if isinstance(obj, h5py.Dataset):
        out["shape"] = [int(x) for x in obj.shape]
        out["dtype"] = str(obj.dtype)
        arr = _sample_numeric_from_dataset(obj)
        out["sample_size"] = int(arr.size)
        out["sample_integer_like_fraction"] = _integer_like_fraction(arr)
        out["sample_min"] = float(np.min(arr)) if arr.size else None
        out["sample_max"] = float(np.max(arr)) if arr.size else None
        return out
    if isinstance(obj, h5py.Group):
        out["attrs"] = {k: (v.tolist() if hasattr(v, "tolist") else str(v)) for k, v in obj.attrs.items()}
        out["keys"] = sorted(list(obj.keys()))
        if "shape" in obj.attrs:
            shp = np.asarray(obj.attrs["shape"])
            out["shape"] = [int(x) for x in shp.tolist()]
        if "data" in obj and isinstance(obj["data"], h5py.Dataset):
            ds = obj["data"]
            out["data_dtype"] = str(ds.dtype)
            arr = _sample_numeric_from_dataset(ds)
            out["sample_size"] = int(arr.size)
            out["sample_integer_like_fraction"] = _integer_like_fraction(arr)
            out["sample_min"] = float(np.min(arr)) if arr.size else None
            out["sample_max"] = float(np.max(arr)) if arr.size else None
        return out
    return out


def _count_source_guess(arch: dict[str, Any]) -> dict[str, Any]:
    layer_info = arch.get("layers", {})
    for name in layer_info:
        lname = str(name).lower()
        if "count" in lname or lname in {"raw", "umi", "umis"}:
            return {"source": f"layers/{name}", "reason": "counts-like layer name"}

    if arch.get("raw_X") is not None:
        frac = arch["raw_X"].get("sample_integer_like_fraction")
        if frac is None or frac >= 0.95:
            return {"source": "raw/X", "reason": "raw slot present and integer-like sample"}

    x = arch.get("X")
    if x is not None:
        frac = x.get("sample_integer_like_fraction")
        if frac is not None and frac >= 0.95:
            return {"source": "X", "reason": "X sample appears integer-like"}

    return {"source": None, "reason": "no clear raw-count source from naming/sample"}


def inspect_architecture(f: h5py.File) -> dict[str, Any]:
    arch: dict[str, Any] = {
        "top_keys": sorted(list(f.keys())),
        "X": None,
        "layers": {},
        "raw_present": bool("raw" in f),
        "raw_X": None,
    }
    if "X" in f:
        arch["X"] = _inspect_matrix_obj(f["X"])
    if "layers" in f and isinstance(f["layers"], h5py.Group):
        for lname in sorted(list(f["layers"].keys())):
            arch["layers"][lname] = _inspect_matrix_obj(f["layers"][lname])
    if "raw" in f and isinstance(f["raw"], h5py.Group) and "X" in f["raw"]:
        arch["raw_X"] = _inspect_matrix_obj(f["raw"]["X"])
    arch["raw_counts_guess"] = _count_source_guess(arch)
    return arch


def inspect_dataset(
    spec: DatasetSpec,
    pert_cols: list[str],
    control_labels: list[str],
) -> dict[str, Any]:
    out: dict[str, Any] = {
        "dataset": spec.name,
        "path": str(spec.h5ad),
        "exists": spec.h5ad.exists(),
        "n_obs": 0,
        "n_vars": 0,
        "obs_columns": [],
        "var_columns": [],
        "candidate_pert_cols": [],
        "recommendation": {
            "pert_col": None,
            "control_label": None,
            "control_count": 0,
            "ready": False,
        },
        "architecture": None,
        "error": None,
    }
    if not spec.h5ad.exists():
        return out

    try:
        with h5py.File(spec.h5ad, "r") as f:
            out["architecture"] = inspect_architecture(f)
            obs_keys = sorted(list(f["obs"].keys())) if "obs" in f else []
            var_keys = sorted(list(f["var"].keys())) if "var" in f else []
            out["obs_columns"] = [k for k in obs_keys if k != "_index"]
            out["var_columns"] = [k for k in var_keys if k != "_index"]
            # AnnData v0.7+: obs/_index dataset exists.
            if "obs" in f and "_index" in f["obs"]:
                out["n_obs"] = int(len(f["obs"]["_index"]))
            # Legacy dataframe encoding: obs attrs contain _index column name.
            elif "obs" in f and "_index" in f["obs"].attrs:
                obs_index_col = str(f["obs"].attrs["_index"])
                if obs_index_col in f["obs"]:
                    out["n_obs"] = int(len(f["obs"][obs_index_col]))

            if "var" in f and "_index" in f["var"]:
                out["n_vars"] = int(len(f["var"]["_index"]))
            elif "var" in f and "_index" in f["var"].attrs:
                var_index_col = str(f["var"].attrs["_index"])
                if var_index_col in f["var"]:
                    out["n_vars"] = int(len(f["var"][var_index_col]))

            candidates: list[dict[str, Any]] = []
            for col in pert_cols:
                c = inspect_obs_column(f, col, control_labels)
                if c.get("present", False):
                    candidates.append(c)
            out["candidate_pert_cols"] = candidates
    except Exception as e:
        out["error"] = f"{type(e).__name__}: {e}"
        return out

    # Pick recommendation by highest control count among candidate labels.
    best_col = None
    best_label = None
    best_count = -1
    for c in out["candidate_pert_cols"]:
        for lbl, cnt in c.get("control_counts", {}).items():
            if int(cnt) > best_count:
                best_count = int(cnt)
                best_col = c["name"]
                best_label = lbl

    if best_col is None and out["candidate_pert_cols"]:
        best_col = out["candidate_pert_cols"][0]["name"]
        best_label = control_labels[0]
        best_count = 0

    out["recommendation"] = {
        "pert_col": best_col,
        "control_label": best_label,
        "control_count": max(best_count, 0),
        "ready": bool(best_col is not None and best_count > 0),
    }
    return out


def load_specs_from_config(config_path: Path) -> list[DatasetSpec]:
    if yaml is None:
        raise RuntimeError("PyYAML is not installed but --config was provided.")
    cfg = yaml.safe_load(config_path.read_text()) or {}
    ds = cfg.get("datasets", {})
    specs: list[DatasetSpec] = []
    for name, dcfg in ds.items():
        if not isinstance(dcfg, dict):
            continue
        if not dcfg.get("enabled", True):
            continue
        h5ad = dcfg.get("h5ad")
        if h5ad:
            specs.append(DatasetSpec(name=name, h5ad=Path(h5ad)))
    return specs


def parse_dataset_args(items: list[str]) -> list[DatasetSpec]:
    specs: list[DatasetSpec] = []
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --dataset '{item}'. Use NAME=path/to/file.h5ad")
        name, path = item.split("=", 1)
        specs.append(DatasetSpec(name=name.strip(), h5ad=Path(path.strip())))
    return specs


def run_analysis(args: argparse.Namespace) -> None:
    specs: list[DatasetSpec] = []
    if args.config is not None:
        specs.extend(load_specs_from_config(args.config))
    if args.use_defaults:
        specs.extend([DatasetSpec(name=k, h5ad=Path(v)) for k, v in DEFAULT_DATASETS.items()])
    if args.dataset:
        specs.extend(parse_dataset_args(args.dataset))

    # De-duplicate by dataset name (last one wins).
    dedup: dict[str, DatasetSpec] = {}
    for s in specs:
        dedup[s.name] = s
    specs = list(dedup.values())
    if not specs:
        raise SystemExit("No datasets provided. Use --config and/or --use-defaults and/or --dataset.")

    pert_cols = CANDIDATE_PERT_COLS + args.candidate_pert_col
    control_labels = CANDIDATE_CONTROL_LABELS + args.candidate_control_label

    args.outdir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    for spec in specs:
        print(f"[inspect] {spec.name}: {spec.h5ad}")
        r = inspect_dataset(spec, pert_cols, control_labels)
        results.append(r)
        (args.outdir / f"{spec.name}.json").write_text(json.dumps(r, indent=2))

    summary_rows = []
    for r in results:
        rec = r.get("recommendation", {})
        arch_guess = (r.get("architecture") or {}).get("raw_counts_guess", {})
        summary_rows.append(
            {
                "dataset": r["dataset"],
                "path": r["path"],
                "exists": bool(r["exists"]),
                "n_obs": int(r.get("n_obs", 0)),
                "n_vars": int(r.get("n_vars", 0)),
                "recommended_pert_col": rec.get("pert_col"),
                "recommended_control_label": rec.get("control_label"),
                "recommended_control_count": int(rec.get("control_count", 0)),
                "ready": bool(rec.get("ready", False)),
                "raw_counts_source_guess": arch_guess.get("source"),
                "raw_counts_reason": arch_guess.get("reason"),
                "error": r.get("error"),
            }
        )

    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.outdir / "summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    (args.outdir / "all_results.json").write_text(json.dumps(results, indent=2))

    print("\n=== Summary ===")
    if not summary_df.empty:
        print(summary_df.to_string(index=False))
    print(f"\nWrote: {summary_path}")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=None, help="Optional config.yaml to load enabled datasets.")
    ap.add_argument(
        "--dataset",
        action="append",
        default=[],
        help="Extra dataset spec NAME=path/to/file.h5ad (can repeat).",
    )
    ap.add_argument(
        "--use-defaults",
        action="store_true",
        help="Include built-in defaults for HCT116/HEK293T/STEMCELL/K562_GWPS_RAW_FULL.",
    )
    ap.add_argument(
        "--candidate-pert-col",
        action="append",
        default=[],
        help="Extra candidate pert_col name(s).",
    )
    ap.add_argument(
        "--candidate-control-label",
        action="append",
        default=[],
        help="Extra candidate control label(s).",
    )
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    dataset = str(snk.params.dataset)
    h5ad = str(snk.input.h5ad)
    return argparse.Namespace(
        config=None,
        dataset=[f"{dataset}={h5ad}"],
        use_defaults=False,
        candidate_pert_col=[],
        candidate_control_label=[],
        outdir=Path(str(snk.params.outdir)),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
