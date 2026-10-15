#!/usr/bin/env python3
"""Build chunk manifests for large perturb-seq datasets.

Supports two sources:
- guide-call CSV with columns: cell_barcode, gene_target
- h5ad with perturbation labels in adata.obs[pert_col]

Each chunk has:
- a subset of perturbations
- a copy (optionally downsampled) of control cells
"""


import argparse
import gzip
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
import sys
sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).resolve().parent))
from workflow_support import require_unique, require_barcodes, read_tables, _open_minimal_backed_adata, _close_backing_resources
from typing import Dict, Iterable, List

import pandas as pd


def split_chunks(items: List[str], size: int) -> List[List[str]]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def split_weighted_chunks(
    items: List[str], counts: Counter, max_items: int, max_cells: int
) -> List[List[str]]:
    """Balance target abundance while respecting target and cell caps."""
    if max_cells <= 0:
        return split_chunks(items, max_items)
    groups: List[List[str]] = []
    loads: List[int] = []
    for item in sorted(items, key=lambda value: (-int(counts[value]), value)):
        weight = int(counts[item])
        candidates = [
            idx for idx, group in enumerate(groups)
            if len(group) < max_items and loads[idx] + weight <= max_cells
        ]
        if candidates:
            idx = min(candidates, key=lambda value: (loads[value], len(groups[value]), value))
            groups[idx].append(item)
            loads[idx] += weight
        else:
            groups.append([item])
            loads.append(weight)
    for group in groups:
        group.sort()
    return groups


def iter_csv_chunks(
    csv_path: Path, chunk_size: int, exclude_guides: set[str] | None = None
) -> Iterable[pd.DataFrame]:
    """Yield guide calls after applying configured guide-level exclusions.

    Exclusions are applied before counting controls or assigning chunks, so an
    excluded control cannot enter a reference set by accident.
    """
    exclude_guides = exclude_guides or set()
    columns = pd.read_csv(csv_path, nrows=0, compression="infer").columns.tolist()
    guide_col = next((name for name in ("guide", "guide_name") if name in columns), None)
    if exclude_guides and guide_col is None:
        raise RuntimeError(
            f"guide exclusions were configured but {csv_path} has neither 'guide' nor 'guide_name'"
        )
    usecols = ["cell_barcode", "gene_target"] + ([guide_col] if guide_col else [])
    with pd.read_csv(
        csv_path,
        usecols=usecols,
        dtype={"cell_barcode": "string", "gene_target": "string"},
        chunksize=chunk_size,
        compression="infer",
        low_memory=True,
    ) as reader:
        for ch in reader:
            if guide_col:
                ch[guide_col] = ch[guide_col].astype("string")
                ch = ch.loc[~ch[guide_col].isin(exclude_guides)].copy()
            yield ch[["cell_barcode", "gene_target"]]


def choose_control_subset(control_barcodes: List[str], max_controls_per_chunk: int, seed: int) -> List[str]:
    control_barcodes = sorted(control_barcodes)
    if max_controls_per_chunk <= 0 or len(control_barcodes) <= max_controls_per_chunk:
        return control_barcodes
    rng = random.Random(seed)
    return rng.sample(control_barcodes, max_controls_per_chunk)


def load_obs_from_h5ad(h5ad_path: Path, pert_col: str) -> pd.DataFrame:
    import h5py
    from anndata.io import read_elem

    with h5py.File(h5ad_path, "r") as handle:
        obs = read_elem(handle["obs"])
    if pert_col not in obs:
        raise ValueError(f"Column {pert_col!r} not found in H5AD obs")
    return pd.DataFrame({"cell_barcode": obs.index.astype(str),
                         "gene_target": obs[pert_col].astype("string").to_numpy()})


def build_from_dataframe(
    df: pd.DataFrame,
    outdir: Path,
    control_label: str,
    perturbations_per_chunk: int,
    max_perturbed_cells_per_chunk: int,
    min_cells_per_perturbation: int,
    max_controls_per_chunk: int,
    control_sample_seed: int,
) -> dict:
    require_unique(df, ["cell_barcode"], "guide calls")
    if df["gene_target"].isna().any() or df["gene_target"].astype(str).str.strip().eq("").any():
        raise ValueError("Guide calls contain missing perturbation labels")
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

    if not control_barcodes:
        raise ValueError(f"No controls found for label {control_label!r}")
    if not perturbations:
        raise ValueError("No perturbations meet min_cells_per_perturbation")
    groups = split_weighted_chunks(
        perturbations, gene_counts, perturbations_per_chunk, max_perturbed_cells_per_chunk
    )
    pert_to_chunk: Dict[str, int] = {}
    for cid, group in enumerate(groups):
        for p in group:
            pert_to_chunk[p] = cid

    controls_for_chunk = choose_control_subset(
        control_barcodes=control_barcodes,
        max_controls_per_chunk=max_controls_per_chunk,
        seed=control_sample_seed,
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
    max_perturbed_cells_per_chunk: int,
    min_cells_per_perturbation: int,
    max_controls_per_chunk: int,
    chunk_size: int,
    control_sample_seed: int,
    exclude_guides: set[str] | None = None,
) -> dict:
    # Keep CSV mode streaming, because this can be very large.
    total_rows = 0
    gene_counts: Counter = Counter()
    control_barcodes: List[str] = []

    seen = set()
    for ch in iter_csv_chunks(csv_path, chunk_size, exclude_guides):
        require_unique(ch, ["cell_barcode"], "guide calls")
        if ch["gene_target"].isna().any() or ch["gene_target"].str.strip().eq("").any():
            raise ValueError("Guide calls contain missing perturbation labels")
        barcodes = set(ch["cell_barcode"].astype(str))
        if seen.intersection(barcodes):
            raise ValueError("Guide calls contain duplicate cell barcodes across CSV blocks")
        seen.update(barcodes)
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
    if not control_barcodes:
        raise ValueError(f"No controls found for label {control_label!r}")
    if not perturbations:
        raise ValueError("No perturbations meet min_cells_per_perturbation")
    groups = split_weighted_chunks(
        perturbations, gene_counts, perturbations_per_chunk, max_perturbed_cells_per_chunk
    )

    pert_to_chunk: Dict[str, int] = {}
    for cid, group in enumerate(groups):
        for p in group:
            pert_to_chunk[p] = cid

    controls_for_chunk = choose_control_subset(
        control_barcodes=control_barcodes,
        max_controls_per_chunk=max_controls_per_chunk,
        seed=control_sample_seed,
    )

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

        for ch in iter_csv_chunks(csv_path, chunk_size, exclude_guides):
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


def run_analysis(args: argparse.Namespace) -> None:
    args.csv = Path(args.csv) if args.csv is not None else None
    args.h5ad = Path(args.h5ad) if args.h5ad is not None else None
    args.outdir = Path(args.outdir)
    args.exclude_guides = {str(x) for x in args.exclude_guides if str(x)}

    if bool(args.csv) == bool(args.h5ad):
        raise RuntimeError("Provide exactly one source: either --csv or --h5ad")

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    if args.exclude_guides and not args.csv:
        raise ValueError("Guide exclusions require guide_calls_csv with guide identifiers")

    if args.csv:
        result = build_from_csv(
            csv_path=args.csv,
            outdir=outdir,
            control_label=args.control_label,
            perturbations_per_chunk=args.perturbations_per_chunk,
            max_perturbed_cells_per_chunk=args.max_perturbed_cells_per_chunk,
            min_cells_per_perturbation=args.min_cells_per_perturbation,
            max_controls_per_chunk=args.max_controls_per_chunk,
            chunk_size=args.chunk_size,
            control_sample_seed=args.control_sample_seed,
            exclude_guides=args.exclude_guides,
        )
        source_meta = {
            "source_type": "csv",
            "csv": str(args.csv),
            "excluded_guides": sorted(args.exclude_guides),
        }
    else:
        obs_df = load_obs_from_h5ad(args.h5ad, args.pert_col)
        result = build_from_dataframe(
            df=obs_df,
            outdir=outdir,
            control_label=args.control_label,
            perturbations_per_chunk=args.perturbations_per_chunk,
            max_perturbed_cells_per_chunk=args.max_perturbed_cells_per_chunk,
            min_cells_per_perturbation=args.min_cells_per_perturbation,
            max_controls_per_chunk=args.max_controls_per_chunk,
            control_sample_seed=args.control_sample_seed,
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
        "max_perturbed_cells_per_chunk": args.max_perturbed_cells_per_chunk,
        "min_cells_per_perturbation": args.min_cells_per_perturbation,
        "max_controls_per_chunk": args.max_controls_per_chunk,
        "control_sample_seed": args.control_sample_seed,
        **result,
    }

    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {outdir / 'manifest.json'}")
    print(f"Wrote {result['chunk_summary_tsv']}")
    print(f"n_chunks={result['n_chunks']}")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path)
    ap.add_argument("--h5ad", type=Path)
    ap.add_argument("--pert-col", type=str, default="gene_target")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--control-label", type=str, default="Non-Targeting")
    ap.add_argument("--perturbations-per-chunk", type=int, default=128)
    ap.add_argument("--max-perturbed-cells-per-chunk", type=int, default=0)
    ap.add_argument("--min-cells-per-perturbation", type=int, default=30)
    ap.add_argument("--max-controls-per-chunk", type=int, default=50000)
    ap.add_argument("--control-sample-seed", type=int, default=0)
    ap.add_argument("--chunk-size", type=int, default=300_000)
    ap.add_argument("--exclude-guides", nargs="*", default=[])
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    source_path = Path(str(snk.input.source))
    source_type = str(snk.params.source_type)
    csv_path = source_path if source_type == "csv" else None
    h5ad_path = source_path if source_type != "csv" else None
    return argparse.Namespace(
        csv=csv_path,
        h5ad=h5ad_path,
        pert_col=str(snk.params.pert_col),
        outdir=Path(str(Path(str(snk.output.manifest)).parent)),
        control_label=str(snk.params.control),
        perturbations_per_chunk=int(snk.params.ppc),
        max_perturbed_cells_per_chunk=int(snk.params.max_perturbed_cells),
        min_cells_per_perturbation=int(snk.params.min_cells),
        max_controls_per_chunk=int(snk.params.max_controls),
        control_sample_seed=int(snk.params.control_seed),
        chunk_size=int(snk.params.read_chunk),
        exclude_guides=list(snk.params.exclude_guides),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
