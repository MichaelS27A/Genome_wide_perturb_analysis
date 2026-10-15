#!/usr/bin/env python3
"""Merge chunk-level global DE outputs into dataset-level tables."""


import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).resolve().parent))
from workflow_support import require_unique, require_barcodes, read_tables, _open_minimal_backed_adata, _close_backing_resources

import pandas as pd


def _read_many(paths: list[Path]) -> list[pd.DataFrame]:
    return read_tables(paths)


def _merge_de_tables(frames: list[pd.DataFrame]) -> pd.DataFrame:
    if not frames:
        return pd.DataFrame(
            columns=[
                "perturbation",
                "rank",
                "gene",
                "score",
                "logfoldchange",
                "pval",
                "pval_adj",
                "pct_nz_group",
                "pct_nz_reference",
            ]
        )
    de = pd.concat(frames, axis=0, ignore_index=True)
    dedupe_cols = [c for c in ("perturbation", "gene", "rank") if c in de.columns]
    if dedupe_cols:
        require_unique(de, dedupe_cols, "merged DE table")
    else:
        de = de.drop_duplicates()
    return de


def _merge_group_counts(frames: list[pd.DataFrame]) -> tuple[pd.DataFrame, str]:
    if not frames:
        return pd.DataFrame(columns=["perturbation", "n_cells"]), "perturbation"

    counts = pd.concat(frames, axis=0, ignore_index=True)
    if "n_cells" not in counts.columns:
        return pd.DataFrame(columns=["perturbation", "n_cells"]), "perturbation"

    group_cols = [c for c in counts.columns if c != "n_cells"]
    group_col = group_cols[0] if group_cols else "perturbation"
    if group_col not in counts.columns:
        counts[group_col] = pd.NA

    counts[group_col] = counts[group_col].astype(str)
    counts["n_cells"] = pd.to_numeric(counts["n_cells"], errors="coerce").fillna(0).astype(int)

    # Controls are repeated across chunks in chunk manifests, so we merge using max
    # rather than sum to avoid artificial inflation.
    merged = (
        counts.groupby(group_col, as_index=False)["n_cells"]
        .max()
        .sort_values("n_cells", ascending=False)
        .reset_index(drop=True)
    )
    return merged, group_col


def _read_chunk_meta(meta_paths: list[Path]) -> list[dict]:
    return [json.loads(p.read_text()) for p in meta_paths]


def run_analysis(args: argparse.Namespace) -> None:
    args.chunk_de = [Path(str(p)) for p in args.chunk_de]
    args.chunk_de_full = [Path(str(p)) for p in args.chunk_de_full]
    args.chunk_group_counts = [Path(str(p)) for p in args.chunk_group_counts]
    args.chunk_meta = [Path(str(p)) for p in args.chunk_meta]
    args.outdir = Path(args.outdir)

    args.outdir.mkdir(parents=True, exist_ok=True)

    de_frames = _read_many(args.chunk_de)
    de_full_frames = _read_many(args.chunk_de_full)
    count_frames = _read_many(args.chunk_group_counts)
    chunk_meta = _read_chunk_meta(args.chunk_meta)

    de = _merge_de_tables(de_frames)
    de_full = _merge_de_tables(de_full_frames)
    group_counts, group_col = _merge_group_counts(count_frames)

    de_out = args.outdir / "perturbation_differential_genes.tsv.gz"
    de_full_out = args.outdir / "perturbation_differential_genes_full.tsv.gz"
    group_counts_out = args.outdir / "group_cell_counts.tsv.gz"
    meta_out = args.outdir / "de_meta.json"
    done_out = args.outdir / "done.txt"

    de.to_csv(de_out, sep="\t", index=False, compression="gzip")
    de_full.to_csv(de_full_out, sep="\t", index=False, compression="gzip")
    group_counts.to_csv(group_counts_out, sep="\t", index=False, compression="gzip")
    done_out.write_text("ok\n")

    meta = {
        "dataset": args.dataset,
        "method": "global_de_chunked_merge",
        "n_chunk_de_files": len(args.chunk_de),
        "n_chunk_de_full_files": len(args.chunk_de_full),
        "n_chunk_group_count_files": len(args.chunk_group_counts),
        "n_chunk_meta_files": len(args.chunk_meta),
        "n_chunk_meta_loaded": len(chunk_meta),
        "n_de_rows": int(de.shape[0]),
        "n_de_full_rows": int(de_full.shape[0]),
        "n_groups_merged": int(group_counts.shape[0]),
        "group_count_key": group_col,
        "chunk_meta": chunk_meta,
        "outputs": {
            "de": str(de_out),
            "de_full": str(de_full_out),
            "group_counts": str(group_counts_out),
            "done": str(done_out),
        },
    }
    meta_out.write_text(json.dumps(meta, indent=2))


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chunk-de", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-de-full", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-group-counts", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-meta", nargs="+", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        dataset=str(snk.params.dataset),
        chunk_de=[Path(str(p)) for p in list(snk.input.de)],
        chunk_de_full=[Path(str(p)) for p in list(snk.input.de_full)],
        chunk_group_counts=[Path(str(p)) for p in list(snk.input.group_counts)],
        chunk_meta=[Path(str(p)) for p in list(snk.input.chunk_meta)],
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
