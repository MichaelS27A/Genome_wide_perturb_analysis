#!/usr/bin/env python3
"""Run ORA on up/down DE genes per perturbation from a DEG table."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    if pvals.size == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    n = float(len(ranked))
    adj = np.empty_like(ranked, dtype=np.float64)
    prev = 1.0
    for i in range(len(ranked) - 1, -1, -1):
        rank = float(i + 1)
        val = (ranked[i] * n) / rank
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty_like(adj, dtype=np.float64)
    out[order] = np.minimum(adj, 1.0)
    return out


def parse_gmt_files(gmt_files: list[Path]) -> dict[str, set[str]]:
    terms: dict[str, set[str]] = {}
    for gmt in gmt_files:
        source = gmt.stem
        with gmt.open("r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                term = parts[0].strip()
                genes = {g.strip() for g in parts[2:] if g.strip()}
                if not genes:
                    continue
                key = f"{source}:{term}"
                if key in terms:
                    terms[key].update(genes)
                else:
                    terms[key] = set(genes)
    return terms


def read_universe_from_h5ad(h5ad_path: Path) -> set[str]:
    # Backed mode keeps memory usage low while reading var names.
    import anndata as ad

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        return {str(x) for x in adata.var_names.astype(str).tolist()}
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()


def run_ora(
    query_genes: set[str],
    terms: dict[str, set[str]],
    universe: set[str],
    min_term_size: int,
    max_term_size: int,
) -> pd.DataFrame:
    U = int(len(universe))
    n = int(len(query_genes))
    if U == 0 or n == 0:
        return pd.DataFrame(
            columns=[
                "term",
                "k_overlap",
                "query_size",
                "term_size",
                "universe_size",
                "expected_overlap",
                "enrichment_ratio",
                "odds_ratio",
                "p_value",
                "hit_genes",
            ]
        )

    rows = []
    for term, genes in terms.items():
        term_genes = genes & universe
        K = int(len(term_genes))
        if K < int(min_term_size) or K > int(max_term_size):
            continue
        hits = sorted(query_genes & term_genes)
        k = int(len(hits))
        if k == 0:
            continue

        # One-sided enrichment p-value.
        p = float(hypergeom.sf(k - 1, U, K, n))
        expected = (n * K) / float(U)
        enrich = (k / expected) if expected > 0 else np.nan

        # Haldane-Anscombe corrected odds ratio for robustness at small counts.
        a = k + 0.5
        b = (n - k) + 0.5
        c = (K - k) + 0.5
        d = (U - K - n + k) + 0.5
        odds_ratio = (a * d) / (b * c) if (b > 0 and c > 0) else np.nan

        rows.append(
            {
                "term": term,
                "k_overlap": k,
                "query_size": n,
                "term_size": K,
                "universe_size": U,
                "expected_overlap": expected,
                "enrichment_ratio": enrich,
                "odds_ratio": odds_ratio,
                "p_value": p,
                "hit_genes": ";".join(hits),
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "term",
                "k_overlap",
                "query_size",
                "term_size",
                "universe_size",
                "expected_overlap",
                "enrichment_ratio",
                "odds_ratio",
                "p_value",
                "hit_genes",
            ]
        )
    out = pd.DataFrame(rows).sort_values("p_value", ascending=True).reset_index(drop=True)
    out["p_adj_bh"] = benjamini_hochberg(pd.to_numeric(out["p_value"], errors="coerce").to_numpy(dtype=np.float64))
    return out


def run_analysis(args: argparse.Namespace) -> None:
    args.de = Path(args.de)
    args.h5ad = Path(args.h5ad)
    args.outdir = Path(args.outdir)
    args.outdir.mkdir(parents=True, exist_ok=True)

    gmt_files = [Path(str(p)) for p in args.gmt_files]
    gmt_files = [p for p in gmt_files if str(p).strip()]
    if not gmt_files:
        raise RuntimeError("No GMT files provided. Set ora.gmt_files in config.")
    missing = [str(p) for p in gmt_files if not p.exists()]
    if missing:
        raise RuntimeError(f"Missing GMT files: {missing}")

    de = pd.read_csv(args.de, sep="\t", compression="infer")
    if de.empty:
        empty = pd.DataFrame(
            columns=[
                "perturbation",
                "direction",
                "term",
                "k_overlap",
                "query_size",
                "term_size",
                "universe_size",
                "expected_overlap",
                "enrichment_ratio",
                "odds_ratio",
                "p_value",
                "p_adj_bh",
                "hit_genes",
            ]
        )
        empty.to_csv(args.outdir / "ora_terms.tsv.gz", sep="\t", index=False, compression="gzip")
        empty.to_csv(args.outdir / "ora_summary.tsv.gz", sep="\t", index=False, compression="gzip")
        meta = {"note": "DE table is empty.", "de": str(args.de)}
        (args.outdir / "ora_meta.json").write_text(json.dumps(meta, indent=2))
        (args.outdir / "done.txt").write_text("ok\n")
        return

    req = {"perturbation", "gene", "logfoldchange"}
    missing_cols = [c for c in req if c not in de.columns]
    if missing_cols:
        raise RuntimeError(f"DE table missing required columns: {missing_cols}")

    pcol = "pval_adj" if "pval_adj" in de.columns else ("pval" if "pval" in de.columns else None)
    if pcol is None:
        raise RuntimeError("DE table must include either 'pval_adj' or 'pval' for ORA filtering.")

    de = de.copy()
    de["gene"] = de["gene"].astype(str)
    de["logfoldchange"] = pd.to_numeric(de["logfoldchange"], errors="coerce")
    de[pcol] = pd.to_numeric(de[pcol], errors="coerce")
    de = de[np.isfinite(de["logfoldchange"]) & np.isfinite(de[pcol])].copy()
    if de.empty:
        raise RuntimeError("No finite logfoldchange/p-value rows found in DE table.")

    universe = read_universe_from_h5ad(args.h5ad)
    terms = parse_gmt_files(gmt_files)
    if not terms:
        raise RuntimeError("No valid terms parsed from GMT files.")

    records = []
    summary_rows = []
    for pert, g in de.groupby("perturbation", sort=True):
        up_genes = set(
            g.loc[
                (g["logfoldchange"] >= float(args.min_abs_logfc)) & (g[pcol] <= float(args.fdr_alpha)),
                "gene",
            ].tolist()
        )
        dn_genes = set(
            g.loc[
                (g["logfoldchange"] <= -float(args.min_abs_logfc)) & (g[pcol] <= float(args.fdr_alpha)),
                "gene",
            ].tolist()
        )

        for direction, qgenes in (("up", up_genes), ("down", dn_genes)):
            qgenes = qgenes & universe
            qsize = len(qgenes)
            if qsize < int(args.min_deg_genes):
                summary_rows.append(
                    {
                        "perturbation": str(pert),
                        "direction": direction,
                        "n_query_genes": int(qsize),
                        "n_terms_tested": 0,
                        "n_terms_fdr_sig": 0,
                        "best_term": np.nan,
                        "best_p_adj_bh": np.nan,
                    }
                )
                continue

            ora = run_ora(
                query_genes=qgenes,
                terms=terms,
                universe=universe,
                min_term_size=int(args.min_term_size),
                max_term_size=int(args.max_term_size),
            )
            if ora.empty:
                summary_rows.append(
                    {
                        "perturbation": str(pert),
                        "direction": direction,
                        "n_query_genes": int(qsize),
                        "n_terms_tested": 0,
                        "n_terms_fdr_sig": 0,
                        "best_term": np.nan,
                        "best_p_adj_bh": np.nan,
                    }
                )
                continue

            ora = ora.sort_values(["p_adj_bh", "p_value"], ascending=[True, True]).reset_index(drop=True)
            n_terms_tested = int(ora.shape[0])
            n_fdr = int((ora["p_adj_bh"] <= float(args.fdr_alpha)).sum())

            keep = ora.head(int(args.max_terms_per_direction)).copy()
            keep.insert(0, "direction", direction)
            keep.insert(0, "perturbation", str(pert))
            records.append(keep)

            summary_rows.append(
                {
                    "perturbation": str(pert),
                    "direction": direction,
                    "n_query_genes": int(qsize),
                    "n_terms_tested": n_terms_tested,
                    "n_terms_fdr_sig": n_fdr,
                    "best_term": str(ora.iloc[0]["term"]),
                    "best_p_adj_bh": float(ora.iloc[0]["p_adj_bh"]),
                }
            )

    if records:
        terms_df = pd.concat(records, axis=0, ignore_index=True)
    else:
        terms_df = pd.DataFrame(
            columns=[
                "perturbation",
                "direction",
                "term",
                "k_overlap",
                "query_size",
                "term_size",
                "universe_size",
                "expected_overlap",
                "enrichment_ratio",
                "odds_ratio",
                "p_value",
                "p_adj_bh",
                "hit_genes",
            ]
        )
    summary_df = pd.DataFrame(summary_rows).sort_values(
        ["n_terms_fdr_sig", "best_p_adj_bh"], ascending=[False, True], na_position="last"
    )

    terms_df.to_csv(args.outdir / "ora_terms.tsv.gz", sep="\t", index=False, compression="gzip")
    summary_df.to_csv(args.outdir / "ora_summary.tsv.gz", sep="\t", index=False, compression="gzip")

    meta = {
        "de": str(args.de),
        "h5ad": str(args.h5ad),
        "gmt_files": [str(p) for p in gmt_files],
        "n_terms_loaded": int(len(terms)),
        "universe_size": int(len(universe)),
        "fdr_alpha": float(args.fdr_alpha),
        "min_abs_logfc": float(args.min_abs_logfc),
        "min_deg_genes": int(args.min_deg_genes),
        "min_term_size": int(args.min_term_size),
        "max_term_size": int(args.max_term_size),
        "max_terms_per_direction": int(args.max_terms_per_direction),
        "n_ora_rows": int(terms_df.shape[0]),
        "n_summary_rows": int(summary_df.shape[0]),
    }
    (args.outdir / "ora_meta.json").write_text(json.dumps(meta, indent=2))
    (args.outdir / "done.txt").write_text("ok\n")


def parse_cli_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--de", type=Path, required=True)
    ap.add_argument("--h5ad", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--gmt-files", nargs="+", default=[])
    ap.add_argument("--fdr-alpha", type=float, default=0.05)
    ap.add_argument("--min-abs-logfc", type=float, default=0.25)
    ap.add_argument("--min-deg-genes", type=int, default=10)
    ap.add_argument("--max-terms-per-direction", type=int, default=50)
    ap.add_argument("--min-term-size", type=int, default=5)
    ap.add_argument("--max-term-size", type=int, default=5000)
    return ap.parse_args()


def args_from_snakemake(snk) -> argparse.Namespace:
    gmt_files = list(getattr(snk.params, "gmt_files", []))
    return argparse.Namespace(
        de=Path(str(snk.input.de)),
        h5ad=Path(str(snk.input.h5ad)),
        outdir=Path(str(snk.params.outdir)),
        gmt_files=[str(x) for x in gmt_files],
        fdr_alpha=float(snk.params.fdr_alpha),
        min_abs_logfc=float(snk.params.min_abs_logfc),
        min_deg_genes=int(snk.params.min_deg_genes),
        max_terms_per_direction=int(snk.params.max_terms_per_direction),
        min_term_size=int(snk.params.min_term_size),
        max_term_size=int(snk.params.max_term_size),
    )


def main() -> None:
    if "snakemake" in globals():
        args = args_from_snakemake(snakemake)
    else:
        args = parse_cli_args()
    run_analysis(args)


if __name__ == "__main__":
    main()
