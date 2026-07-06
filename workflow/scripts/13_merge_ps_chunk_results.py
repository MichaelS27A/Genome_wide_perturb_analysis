#!/usr/bin/env python3
"""Merge chunk-level PS outputs into dataset-level tables."""


import argparse
import json
from pathlib import Path

import pandas as pd


def _read_many(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for p in paths:
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p, sep="\t", compression="infer")
        except pd.errors.EmptyDataError:
            continue
        if df is None or df.empty:
            continue
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, axis=0, ignore_index=True)


def run_analysis(args: argparse.Namespace) -> None:
    args.chunk_cell_scores = [Path(str(p)) for p in args.chunk_cell_scores]
    args.chunk_summary = [Path(str(p)) for p in args.chunk_summary]
    args.outdir = Path(args.outdir)

    args.outdir.mkdir(parents=True, exist_ok=True)

    cell = _read_many(args.chunk_cell_scores)
    summ = _read_many(args.chunk_summary)

    if not cell.empty and {"cell_barcode", "perturbation"}.issubset(cell.columns):
        cell = cell.drop_duplicates(subset=["cell_barcode", "perturbation"], keep="first")

    if not summ.empty:
        if "perturbation" in summ.columns:
            summ = summ.drop_duplicates(subset=["perturbation"], keep="first")
        else:
            summ = summ.drop_duplicates()

    cell_out = args.outdir / "cell_scores.tsv.gz"
    summ_out = args.outdir / "perturbation_summary.tsv.gz"
    meta_out = args.outdir / "method_meta.json"
    done_out = args.outdir / "done.txt"

    cell.to_csv(cell_out, sep="\t", index=False, compression="gzip")
    summ.to_csv(summ_out, sep="\t", index=False, compression="gzip")

    meta = {
        "dataset": args.dataset,
        "method": "PS_scMAGeCK_chunked",
        "n_chunk_cell_score_files": len(args.chunk_cell_scores),
        "n_chunk_summary_files": len(args.chunk_summary),
        "n_cell_rows": int(cell.shape[0]),
        "n_summary_rows": int(summ.shape[0]),
        "outputs": {
            "cell_scores": str(cell_out),
            "perturbation_summary": str(summ_out),
        },
    }
    meta_out.write_text(json.dumps(meta, indent=2))
    done_out.write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chunk-cell-scores", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-summary", nargs="+", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        dataset=str(snk.params.dataset),
        chunk_cell_scores=[Path(str(p)) for p in list(snk.input.cell_scores)],
        chunk_summary=[Path(str(p)) for p in list(snk.input.summary)],
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
