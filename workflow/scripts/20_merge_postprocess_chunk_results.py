#!/usr/bin/env python3
"""Merge chunk-level postprocess outputs into dataset-level tables."""


import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).resolve().parent))
from workflow_support import require_unique, require_barcodes, read_tables, _open_minimal_backed_adata, _close_backing_resources

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


def _read_many(paths: list[Path]) -> list[pd.DataFrame]:
    return read_tables(paths)


def _merge_frames(frames: list[pd.DataFrame], dedupe_cols: list[str] | None = None) -> pd.DataFrame:
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, axis=0, ignore_index=True)
    if dedupe_cols:
        cols = [c for c in dedupe_cols if c in out.columns]
        if cols:
            require_unique(out, cols, "merged postprocess table")
        else:
            out = out.drop_duplicates()
    else:
        out = out.drop_duplicates()
    return out.reset_index(drop=True)


def _merge_pseudobulk(paths: list[Path]) -> ad.AnnData:
    adatas: list[ad.AnnData] = []
    for p in paths:
        cur = ad.read_h5ad(p)
        if cur.n_obs == 0:
            continue
        adatas.append(cur)

    if not adatas:
        return ad.AnnData()

    if any(not cur.var_names.equals(adatas[0].var_names) for cur in adatas):
        raise ValueError("Chunk pseudobulks have different gene identifiers or order")
    merged = ad.concat(adatas, axis=0, join="inner", merge="same", index_unique=None)
    if not merged.obs_names.is_unique:
        raise ValueError("Duplicate perturbations across chunk pseudobulks")
    return merged


def _read_chunk_meta(paths: list[Path]) -> list[dict]:
    return [json.loads(p.read_text()) for p in paths]


def _embed_merged_pseudobulk(pb: ad.AnnData, seed: int = 0) -> pd.DataFrame:
    """Create the sole comparable perturbation UMAP after all chunks merge."""
    if pb.n_obs == 0:
        return pd.DataFrame(columns=["perturbation", "n_selected_cells", "leiden", "umap_1", "umap_2"])

    embedded = pb.copy()
    if min(embedded.n_obs, embedded.n_vars) >= 3:
        sc.pp.highly_variable_genes(embedded, n_top_genes=min(2000, embedded.n_vars), subset=False)
        sc.pp.pca(embedded, n_comps=min(50, embedded.n_obs - 1, embedded.n_vars - 1), mask_var=None, random_state=seed)
        sc.pp.neighbors(embedded, n_neighbors=min(15, embedded.n_obs - 1), random_state=seed)
        sc.tl.umap(embedded, random_state=seed)
        sc.tl.leiden(embedded, key_added="leiden", flavor="igraph", n_iterations=2, random_state=seed)
        coords = embedded.obsm["X_umap"]
    else:
        embedded.obs["leiden"] = pd.NA
        coords = np.full((embedded.n_obs, 2), np.nan)

    out = embedded.obs.copy().rename_axis("perturbation").reset_index()
    out["umap_1"] = coords[:, 0]
    out["umap_2"] = coords[:, 1]
    return out


def run_analysis(args: argparse.Namespace) -> None:
    args.chunk_pseudobulk = [Path(str(p)) for p in args.chunk_pseudobulk]
    args.chunk_umap = [Path(str(p)) for p in args.chunk_umap]
    args.chunk_de = [Path(str(p)) for p in args.chunk_de]
    args.chunk_de_full = [Path(str(p)) for p in args.chunk_de_full]
    args.chunk_long = [Path(str(p)) for p in args.chunk_long]
    args.chunk_gene_long = [Path(str(p)) for p in args.chunk_gene_long]
    args.chunk_meta = [Path(str(p)) for p in args.chunk_meta]
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    pb = _merge_pseudobulk(args.chunk_pseudobulk)
    # Chunk UMAP files are retained as lightweight workflow sentinels, but
    # their coordinates are intentionally not used: chunks have incompatible
    # coordinate systems.  Embed the merged pseudobulk once instead.
    umap = _embed_merged_pseudobulk(pb, args.random_seed)
    de = _merge_frames(_read_many(args.chunk_de), dedupe_cols=["perturbation", "gene", "rank"])
    de_full = _merge_frames(_read_many(args.chunk_de_full), dedupe_cols=["perturbation", "gene", "rank"])
    long_table = _merge_frames(_read_many(args.chunk_long), dedupe_cols=["perturbation"])
    gene_long = _merge_frames(_read_many(args.chunk_gene_long), dedupe_cols=["perturbation", "gene", "rank"])
    for frame in (long_table, gene_long):
        stale = [c for c in ("leiden", "umap_1", "umap_2") if c in frame.columns]
        if stale:
            frame.drop(columns=stale, inplace=True)
    if not long_table.empty:
        long_table = long_table.drop(columns=["n_selected_cells"], errors="ignore").merge(
            umap, on="perturbation", how="left", validate="one_to_one"
        )
    if not gene_long.empty:
        gene_long = gene_long.merge(
            umap[["perturbation", "leiden", "umap_1", "umap_2", "n_selected_cells"]],
            on="perturbation", how="left", validate="many_to_one"
        )
    chunk_meta = _read_chunk_meta(args.chunk_meta)

    pb_out = args.outdir / "perturbation_pseudobulk.h5ad"
    umap_out = args.outdir / "perturbation_umap_leiden.tsv.gz"
    de_out = args.outdir / "perturbation_differential_genes.tsv.gz"
    de_full_out = args.outdir / "perturbation_differential_genes_full.tsv.gz"
    long_out = args.outdir / "perturbation_long_table.tsv.gz"
    gene_long_out = args.outdir / "perturbation_gene_long_table.tsv.gz"
    meta_out = args.outdir / "postprocess_meta.json"

    pb.write_h5ad(pb_out, compression="gzip")
    umap.to_csv(umap_out, sep="\t", index=False, compression="gzip")
    de.to_csv(de_out, sep="\t", index=False, compression="gzip")
    de_full.to_csv(de_full_out, sep="\t", index=False, compression="gzip")
    long_table.to_csv(long_out, sep="\t", index=False, compression="gzip")
    gene_long.to_csv(gene_long_out, sep="\t", index=False, compression="gzip")

    meta = {
        "dataset": args.dataset,
        "method": "postprocess_chunked_merge_with_dataset_level_umap",
        "n_chunk_pseudobulk_files": len(args.chunk_pseudobulk),
        "n_chunk_umap_files": len(args.chunk_umap),
        "n_chunk_de_files": len(args.chunk_de),
        "n_chunk_de_full_files": len(args.chunk_de_full),
        "n_chunk_long_files": len(args.chunk_long),
        "n_chunk_gene_long_files": len(args.chunk_gene_long),
        "n_chunk_meta_files": len(args.chunk_meta),
        "n_chunk_meta_loaded": len(chunk_meta),
        "n_pseudobulk_rows": int(pb.n_obs),
        "n_de_rows": int(de.shape[0]),
        "n_de_full_rows": int(de_full.shape[0]),
        "n_long_rows": int(long_table.shape[0]),
        "n_gene_long_rows": int(gene_long.shape[0]),
        "chunk_meta": chunk_meta,
        "outputs": {
            "pseudobulk": str(pb_out),
            "umap_leiden": str(umap_out),
            "de": str(de_out),
            "de_full": str(de_full_out),
            "long_table": str(long_out),
            "gene_long_table": str(gene_long_out),
        },
    }
    meta_out.write_text(json.dumps(meta, indent=2))


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--chunk-pseudobulk", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-umap", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-de", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-de-full", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-long", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-gene-long", nargs="+", type=Path, required=True)
    ap.add_argument("--chunk-meta", nargs="+", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--random-seed", type=int, default=0)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    return argparse.Namespace(
        dataset=str(snk.params.dataset),
        random_seed=int(snk.params.seed),
        chunk_pseudobulk=[Path(str(p)) for p in list(snk.input.pseudobulk)],
        chunk_umap=[Path(str(p)) for p in list(snk.input.umap_leiden)],
        chunk_de=[Path(str(p)) for p in list(snk.input.de)],
        chunk_de_full=[Path(str(p)) for p in list(snk.input.de_full)],
        chunk_long=[Path(str(p)) for p in list(snk.input.long_table)],
        chunk_gene_long=[Path(str(p)) for p in list(snk.input.gene_long_table)],
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
