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

    if args.require_ready and not ready:
        raise AssertionError(
            f"Dataset '{args.dataset}' not ready from audit: "
            f"pert_col={pert_col!r}, control_label={ctrl!r}, ready={ready}."
        )

    if args.require_raw_source and not raw_guess:
        raise AssertionError(
            f"Dataset '{args.dataset}' has no inferred raw-count source from audit architecture."
        )

    print(
        f"[audit-ok] dataset={args.dataset} ready={ready} pert_col={pert_col} "
        f"control_label={ctrl} raw_source={raw_guess}"
    )
    if args.ready_path is not None:
        args.ready_path.parent.mkdir(parents=True, exist_ok=True)
        args.ready_path.write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit-json", type=Path, required=True)
    ap.add_argument("--dataset", type=str, required=True)
    ap.add_argument("--require-ready", action="store_true")
    ap.add_argument("--require-raw-source", action="store_true")
    ap.add_argument("--ready-path", type=Path, default=None)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        audit_json=Path(str(snk.input.json)),
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
