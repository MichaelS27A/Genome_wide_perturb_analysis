#!/usr/bin/env python3
"""Build chunk manifests for large perturb-seq datasets.

Supports two sources:
- guide-call CSV with columns: cell_barcode, gene_target
- h5ad with perturbation labels in adata.obs[pert_col]

Each chunk has:
- a subset of perturbations
- a copy (optionally downsampled) of control cells
"""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd


def split_chunks(items: List[str], size: int) -> List[List[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def iter_csv_chunks(csv_path: Path, chunk_size: int) -> Iterable[pd.DataFrame]:
    for ch in pd.read_csv(
        csv_path,
        usecols=["cell_barcode", "gene_target"],
        dtype={"cell_barcode": "string", "gene_target": "string"},
        chunksize=chunk_size,
        compression="infer",
        low_memory=True,
    ):
        yield ch


def load_obs_from_h5ad(h5ad_path: Path, pert_col: str) -> pd.DataFrame:
    import scanpy as sc

    adata = sc.read_h5ad(h5ad_path, backed="r")
    try:
        if pert_col not in adata.obs.columns:
            raise RuntimeError(f"Column '{pert_col}' not found in {h5ad_path}.obs")

        obs = adata.obs[[pert_col]].copy()
        obs = obs.reset_index().rename(columns={"index": "cell_barcode", pert_col: "gene_target"})
        obs["cell_barcode"] = obs["cell_barcode"].astype("string")
        obs["gene_target"] = obs["gene_target"].astype("string")
        return obs
    finally:
        adata.file.close()


def build_from_dataframe(
    df: pd.DataFrame,
    outdir: Path,
    control_label: str,
    perturbations_per_chunk: int,
    min_cells_per_perturbation: int,
    max_controls_per_chunk: int,
) -> dict:
    chunks_dir = outdir / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    total_rows = int(df.shape[0])
    gt = df["gene_target"].fillna("<NA>").astype(str)
    gene_counts = Counter(gt.tolist())

    control_barcodes = (
        df.loc[gt == control_label, "cell_barcode"].dropna().astype(str).tolist()
    )

    perturbations = [
        p
        for p, c in gene_counts.items()
        if p not in {control_label, "<NA>"} and c >= min_cells_per_perturbation
    ]
    perturbations.sort()

    groups = split_chunks(perturbations, perturbations_per_chunk)
    pert_to_chunk: Dict[str, int] = {}
    for cid, group in enumerate(groups):
        for p in group:
            pert_to_chunk[p] = cid

    controls_for_chunk = (
        control_barcodes[:max_controls_per_chunk] if max_controls_per_chunk > 0 else control_barcodes
    )

    writers = {}
    summary = defaultdict(lambda: {"n_perturbed_cells": 0, "n_control_cells": 0})

    try:
        for cid, group in enumerate(groups):
            cell_path = chunks_dir / f"chunk_{cid:04d}_cells.tsv.gz"
            w = gzip.open(cell_path, "wt")
            w.write("cell_barcode\tgene_target\tis_control\n")

            for cb in controls_for_chunk:
                w.write(f"{cb}\t{control_label}\t1\n")
            summary[cid]["n_control_cells"] = len(controls_for_chunk)

            pert_path = chunks_dir / f"chunk_{cid:04d}_perturbations.txt"
            pert_path.write_text("\n".join(group) + "\n")
            writers[cid] = w

        df2 = df.dropna(subset=["cell_barcode", "gene_target"]).copy()
        df2["gene_target"] = df2["gene_target"].astype(str)
        df2 = df2[df2["gene_target"].isin(pert_to_chunk)]

        for cb, gt2 in zip(df2["cell_barcode"].astype(str).tolist(), df2["gene_target"].tolist()):
            cid = pert_to_chunk[gt2]
            writers[cid].write(f"{cb}\t{gt2}\t0\n")
            summary[cid]["n_perturbed_cells"] += 1
    finally:
        for w in writers.values():
            w.close()

    summary_tsv = outdir / "chunk_summary.tsv"
    with summary_tsv.open("w") as fh:
        fh.write(
            "chunk_id\tn_perturbations\tn_perturbed_cells\tn_control_cells\tn_total_cells\tfirst_perturbations\n"
        )
        for cid, group in enumerate(groups):
            n_pc = summary[cid]["n_perturbed_cells"]
            n_cc = summary[cid]["n_control_cells"]
            head = ",".join(group[:5]) + ("..." if len(group) > 5 else "")
            fh.write(f"{cid}\t{len(group)}\t{n_pc}\t{n_cc}\t{n_pc+n_cc}\t{head}\n")

    return {
        "total_rows": total_rows,
        "n_control_cells_global": len(control_barcodes),
        "n_perturbations_all": len(gene_counts),
        "n_perturbations_kept": len(perturbations),
        "n_chunks": len(groups),
        "chunk_summary_tsv": str(summary_tsv),
    }


def build_from_csv(
    csv_path: Path,
    outdir: Path,
    control_label: str,
    perturbations_per_chunk: int,
    min_cells_per_perturbation: int,
    max_controls_per_chunk: int,
    chunk_size: int,
) -> dict:
    # Keep CSV mode streaming, because this can be very large.
    total_rows = 0
    gene_counts: Counter = Counter()
    control_barcodes: List[str] = []

    for ch in iter_csv_chunks(csv_path, chunk_size):
        total_rows += len(ch)
        gt = ch["gene_target"].fillna("<NA>")
        vc = gt.value_counts(dropna=False)
        for k, v in vc.items():
            gene_counts[str(k)] += int(v)

        ctl = ch.loc[gt == control_label, "cell_barcode"].dropna().astype(str)
        control_barcodes.extend(ctl.tolist())

    perturbations = [
        p
        for p, c in gene_counts.items()
        if p not in {control_label, "<NA>"} and c >= min_cells_per_perturbation
    ]
    perturbations.sort()
    groups = split_chunks(perturbations, perturbations_per_chunk)

    pert_to_chunk: Dict[str, int] = {}
    for cid, group in enumerate(groups):
        for p in group:
            pert_to_chunk[p] = cid

    controls_for_chunk = control_barcodes[: max_controls_per_chunk] if max_controls_per_chunk > 0 else control_barcodes

    chunks_dir = outdir / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    writers = {}
    summary = defaultdict(lambda: {"n_perturbed_cells": 0, "n_control_cells": 0})

    try:
        for cid, group in enumerate(groups):
            cell_path = chunks_dir / f"chunk_{cid:04d}_cells.tsv.gz"
            w = gzip.open(cell_path, "wt")
            w.write("cell_barcode\tgene_target\tis_control\n")

            for cb in controls_for_chunk:
                w.write(f"{cb}\t{control_label}\t1\n")
            summary[cid]["n_control_cells"] = len(controls_for_chunk)

            pert_path = chunks_dir / f"chunk_{cid:04d}_perturbations.txt"
            pert_path.write_text("\n".join(group) + "\n")
            writers[cid] = w

        for ch in iter_csv_chunks(csv_path, chunk_size):
            ch = ch.dropna(subset=["cell_barcode", "gene_target"])
            ch["gene_target"] = ch["gene_target"].astype(str)
            ch = ch[ch["gene_target"].isin(pert_to_chunk)]

            for cb, gt2 in zip(ch["cell_barcode"].astype(str).tolist(), ch["gene_target"].tolist()):
                cid = pert_to_chunk[gt2]
                writers[cid].write(f"{cb}\t{gt2}\t0\n")
                summary[cid]["n_perturbed_cells"] += 1
    finally:
        for w in writers.values():
            w.close()

    summary_tsv = outdir / "chunk_summary.tsv"
    with summary_tsv.open("w") as fh:
        fh.write(
            "chunk_id\tn_perturbations\tn_perturbed_cells\tn_control_cells\tn_total_cells\tfirst_perturbations\n"
        )
        for cid, group in enumerate(groups):
            n_pc = summary[cid]["n_perturbed_cells"]
            n_cc = summary[cid]["n_control_cells"]
            head = ",".join(group[:5]) + ("..." if len(group) > 5 else "")
            fh.write(f"{cid}\t{len(group)}\t{n_pc}\t{n_cc}\t{n_pc+n_cc}\t{head}\n")

    return {
        "total_rows": total_rows,
        "n_control_cells_global": len(control_barcodes),
        "n_perturbations_all": len(gene_counts),
        "n_perturbations_kept": len(perturbations),
        "n_chunks": len(groups),
        "chunk_summary_tsv": str(summary_tsv),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--h5ad", type=Path)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--perturbations-per-chunk", type=int, default=128)
    ap.add_argument("--min-cells-per-perturbation", type=int, default=30)
    ap.add_argument("--max-controls-per-chunk", type=int, default=50000)
    ap.add_argument("--chunk-size", type=int, default=300_000)
    args = ap.parse_args()

    if bool(args.csv) == bool(args.h5ad):
        raise RuntimeError("Provide exactly one source: either --csv or --h5ad")

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.csv:
        result = build_from_csv(
            csv_path=args.csv,
            outdir=outdir,
            control_label=args.control_label,
            perturbations_per_chunk=args.perturbations_per_chunk,
            min_cells_per_perturbation=args.min_cells_per_perturbation,
            max_controls_per_chunk=args.max_controls_per_chunk,
            chunk_size=args.chunk_size,
        )
        source_meta = {"source_type": "csv", "csv": str(args.csv)}
    else:
        obs_df = load_obs_from_h5ad(args.h5ad, args.pert_col)
        result = build_from_dataframe(
            df=obs_df,
            outdir=outdir,
            control_label=args.control_label,
            perturbations_per_chunk=args.perturbations_per_chunk,
            min_cells_per_perturbation=args.min_cells_per_perturbation,
            max_controls_per_chunk=args.max_controls_per_chunk,
        )
        source_meta = {
            "source_type": "h5ad_obs",
            "h5ad": str(args.h5ad),
            "pert_col": args.pert_col,
        }

    manifest = {
        **source_meta,
        "controls_label": args.control_label,
        "perturbations_per_chunk": args.perturbations_per_chunk,
        "min_cells_per_perturbation": args.min_cells_per_perturbation,
        "max_controls_per_chunk": args.max_controls_per_chunk,
        **result,
    }

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {outdir / 'manifest.json'}")
    print(f"Wrote {result['chunk_summary_tsv']}")
    print(f"n_chunks={result['n_chunks']}")


if __name__ == "__main__":
    main()
