#!/usr/bin/env python3
"""Validate dataset audit JSON and fail fast when dataset is not ready."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit-json", type=Path, required=True)
    ap.add_argument("--dataset", type=str, required=True)
    ap.add_argument("--require-ready", action="store_true")
    ap.add_argument("--require-raw-source", action="store_true")
    args = ap.parse_args()

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


if __name__ == "__main__":
    main()
