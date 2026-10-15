#!/usr/bin/env python3
"""Validate dataset audit JSON and fail fast when dataset is not ready."""


import argparse
import json
from pathlib import Path


def run_analysis(args: argparse.Namespace) -> None:
    args.audit_json = Path(args.audit_json)
    args.ready_path = Path(args.ready_path) if args.ready_path is not None else None

    payload = json.loads(args.audit_json.read_text())
    rec = payload.get("recommendation", {})
    ready = bool(rec.get("ready", False))
    pert_col = rec.get("pert_col")
    ctrl = rec.get("control_label")

    arch = payload.get("architecture") or {}
    raw_guess = (arch.get("raw_counts_guess") or {}).get("source")

    # Validate the configured labels, not the inspector's guessed names.
    if args.h5ad is not None:
        import h5py
        import numpy as np
        import pandas as pd
        from anndata.io import read_elem
        with h5py.File(args.h5ad, "r") as handle:
            obs = read_elem(handle["obs"])
            if not obs.index.is_unique:
                raise ValueError("H5AD observation names must be unique")
            x = handle["X"]
            values = x["data"][:100000] if isinstance(x, h5py.Group) else x[:min(100, x.shape[0]), :min(1000, x.shape[1])]
            if args.require_raw_source and (not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all()):
                raise ValueError("H5AD X must contain raw counts; a counts layer alone is insufficient")
        if args.guide_calls is not None:
            columns = pd.read_csv(args.guide_calls, nrows=0).columns
            if not {"cell_barcode", "gene_target"}.issubset(columns):
                raise ValueError("guide_calls_csv requires cell_barcode,gene_target columns")
            # Full label/barcode validation is performed by the streaming manifest builder.
            ready = True
        else:
            if args.pert_col not in obs or args.control_label not in set(obs[args.pert_col].astype(str)):
                raise ValueError("Configured perturbation column or control label missing from H5AD obs")
            ready = True
    elif args.require_ready and not ready:
        raise ValueError(f"Dataset {args.dataset!r} not ready from audit")
    elif args.require_raw_source and not raw_guess:
        raise ValueError("No inferred raw-count source in audit")

    print(
        f"[audit-ok] dataset={args.dataset} ready={ready} pert_col={pert_col} "
        f"control_label={ctrl} raw_source={raw_guess}"
    )
    if args.ready_path is not None:
        args.ready_path.parent.mkdir(parents=True, exist_ok=True)
        args.ready_path.write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", type=Path)
    ap.add_argument("--guide-calls", type=Path)
    ap.add_argument("--pert-col", default="gene_target")
    ap.add_argument("--control-label", default="Non-Targeting")
    ap.add_argument("--audit-json", type=Path, required=True)
    ap.add_argument("--dataset", type=str, required=True)
    ap.add_argument("--require-ready", action="store_true")
    ap.add_argument("--require-raw-source", action="store_true")
    ap.add_argument("--ready-path", type=Path, default=None)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        audit_json=Path(str(snk.input.json)),
        h5ad=Path(str(snk.input.h5ad)),
        guide_calls=Path(str(snk.input.calls[0])) if snk.input.calls else None,
        pert_col=str(snk.params.pert_col),
        control_label=str(snk.params.control),
        dataset=str(snk.params.dataset),
        require_ready=bool(snk.params.require_ready),
        require_raw_source=bool(snk.params.require_raw),
        ready_path=Path(str(snk.output.ready)),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
