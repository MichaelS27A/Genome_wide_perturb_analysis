"""Optional DESeq2 Wald contrasts on raw pseudobulks with explicit experimental units."""

import importlib.metadata
import json
from pathlib import Path
import sys
import warnings

sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).parent))
from statistics_support import integer_counts, null_guide_sets, stable_seed, versions, write_json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


RESULT_COLUMNS = ["context", "stratum", "perturbation", "guide", "analysis_level", "gene",
                  "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]


def collapse_biological_samples(counts, metadata):
    keys = ["context", "target", "sample"]
    grouped = metadata.groupby(keys, sort=True, observed=True)
    codes = grouped.ngroup().to_numpy()
    assignment = sparse.csr_matrix(
        (np.ones(len(metadata), dtype=np.int64), (codes, np.arange(len(metadata)))),
        shape=(grouped.ngroups, len(metadata)),
    )
    samples = grouped["n_cells"].sum().reset_index()
    samples.index = [f"S{i:08d}" for i in range(len(samples))]
    return (assignment @ counts).tocsr(), samples


def contrast_data(counts, metadata, treatment, reference, mode, paired, settings):
    """Expected ineligible comparisons return a reason; malformed input raises."""
    keep = np.flatnonzero(np.asarray(treatment) | np.asarray(reference))
    if not len(keep):
        return None, None, "no_samples"
    samples = metadata.iloc[keep].copy()
    samples["target"] = np.where(np.asarray(treatment)[keep], "treatment", "control")
    matrix = counts[keep].tocsr()
    if mode == "biological_sample":
        matrix, samples = collapse_biological_samples(matrix, samples)
    eligible = ((samples["n_cells"] >= settings["min_cells_per_pool"])
                & (np.asarray(matrix.sum(axis=1)).ravel() > 0)).to_numpy()
    matrix, samples = matrix[eligible], samples.loc[eligible].copy()
    if mode == "guide_batch":
        samples["block"] = [json.dumps([sample, batch]) for sample, batch in
                            zip(samples["sample"], samples["batch"])]
    block = "sample" if mode == "biological_sample" else "block"
    if paired and mode != "kolf_resample":
        group_blocks = [set(samples.loc[samples["target"].eq(label), block])
                        for label in ("treatment", "control")]
        shared = group_blocks[0] & group_blocks[1]
        keep = samples[block].isin(shared).to_numpy()
        matrix, samples = matrix[keep], samples.loc[keep].copy()
    numbers = samples["target"].value_counts()
    if min(numbers.get("treatment", 0), numbers.get("control", 0)) < settings["min_replicates"]:
        return None, None, "insufficient_matched_samples_or_pools"
    samples.index = [f"S{i:08d}" for i in range(len(samples))]
    return matrix, samples, ""


def fit_wald(counts, samples, genes, settings, block=None):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    keep = ((np.asarray(counts.sum(axis=0)).ravel() >= settings["min_gene_total_count"])
            & (np.asarray((counts > 0).sum(axis=0)).ravel() >= settings["min_gene_pools"]))
    if not keep.any():
        return None, "no_genes_after_count_filter", {}
    counts = counts[:, keep]
    if counts.shape[0] * counts.shape[1] > settings["max_model_elements"]:
        raise ValueError("DESeq2 model exceeds the declared dense element budget")
    design = pd.DataFrame({"Intercept": np.ones(len(samples))}, index=samples.index)
    if block is not None:
        dummies = pd.get_dummies(samples[block], drop_first=True, dtype=float)
        dummies.columns = [f"block_{i}" for i in range(len(dummies.columns))]
        design = pd.concat([design, dummies], axis=1)
    design["treatment"] = samples["target"].eq("treatment").astype(float)
    if np.linalg.matrix_rank(design.to_numpy()) != design.shape[1]:
        return None, "confounded_design", {}
    if len(samples) <= design.shape[1]:
        return None, "no_residual_degrees_of_freedom", {}
    frame = pd.DataFrame(integer_counts(counts).toarray(), index=samples.index,
                         columns=np.asarray(genes)[keep])
    dds = DeseqDataSet(counts=frame, metadata=samples, design=design,
                       size_factors_fit_type="poscounts", fit_type=settings["dispersion_fit"],
                       refit_cooks=True, n_cpus=1, quiet=True)
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message="The dispersion trend curve fitting did not converge.*",
                                category=UserWarning)
        dds.deseq2()
    if dds.uns["disp_function_type"] != settings["dispersion_fit"]:
        raise RuntimeError("DESeq2 could not fit the declared dispersion trend. "
                           "Select a different dispersion_fit explicitly after reviewing the data.")
    contrast = np.zeros(design.shape[1])
    contrast[-1] = 1
    stats = DeseqStats(dds, contrast=contrast, alpha=settings["fdr_alpha"],
                      cooks_filter=True, independent_filter=True, n_cpus=1, quiet=True)
    stats.summary()
    result = stats.results_df.rename_axis("gene").reset_index()
    if not np.isfinite(result["pvalue"]).any():
        raise RuntimeError("DESeq2 produced no finite p-values; inspect the input counts and model")
    return result, "", {
        "n_genes": len(result), "design_columns": design.columns.tolist(),
        "design_rank": design.shape[1],
        "n_finite_pvalues": int(np.isfinite(result["pvalue"]).sum()),
        "n_finite_adjusted_pvalues": int(np.isfinite(result["padj"]).sum()),
        "n_significant": int(result["padj"].lt(settings["fdr_alpha"]).sum()),
        "dispersion_fit": str(dds.uns.get("disp_function_type", "not_recorded")),
    }


def guide_agreement(results, target_genes):
    records = []
    tables = results.loc[results["analysis_level"].eq("guide")]
    for (context, target), frame in tables.groupby(["context", "perturbation"]):
        mapped = target in target_genes
        excluded = target_genes.get(target, set())
        values = frame.loc[~frame["gene"].isin(excluded)].pivot(
            index="gene", columns="guide", values="log2FoldChange")
        for i, first in enumerate(values.columns):
            for second in values.columns[i + 1:]:
                paired = values[[first, second]].replace([np.inf, -np.inf], np.nan).dropna()
                reason = ("unmapped_target_gene" if not mapped else
                          "insufficient_genes" if len(paired) < 3 else "")
                if not reason and (paired.nunique() < 2).any():
                    reason = "constant_effect_vector"
                corr = paired.corr().iloc[0, 1] if not reason else np.nan
                records.append(dict(context=context, perturbation=target, guide_a=first,
                                    guide_b=second, n_genes=len(paired), pearson_r=corr,
                                    target_gene_mapped=mapped, reason=reason))
    return pd.DataFrame(records, columns=["context", "perturbation", "guide_a", "guide_b",
                                         "n_genes", "pearson_r", "target_gene_mapped", "reason"])


def run(counts_path, dataset, settings, outdir):
    source = ad.read_h5ad(counts_path)
    integer_counts(source.X)
    counts, pools = source.X.tocsr(), source.obs.copy()
    symbol_col = dataset["statistics"].get("gene_symbol_col")
    if symbol_col:
        if symbol_col not in source.var:
            raise ValueError(f"Declared gene_symbol_col is missing from var: {symbol_col}")
        target_genes = {}
        for gene, symbol in source.var[symbol_col].dropna().items():
            target_genes.setdefault(str(symbol), set()).add(str(gene))
    else:
        target_genes = {str(gene): {str(gene)} for gene in source.var_names}
    mode = dataset["statistics"]["pseudobulk_mode"]
    paired = dataset["statistics"].get("paired", True) and mode != "kolf_resample"
    control = dataset["control_label"]
    results, summaries = [], []

    def evaluate(frame, indices, treat, reference, target, level, guide=""):
        matrix, samples, reason = contrast_data(
            counts[indices], frame, treat, reference, mode, paired, settings)
        stratum = (json.dumps([frame["sample"].iloc[0], frame["batch"].iloc[0]])
                   if mode == "kolf_resample" else "pooled_samples")
        record = dict(context=frame["context"].iloc[0], stratum=stratum, perturbation=target,
                      analysis_level=level, guide=guide, mode=mode)
        if not reason:
            block = ("sample" if mode == "biological_sample" else "block") if paired else None
            result, reason, diagnostics = fit_wald(matrix, samples, source.var_names, settings, block)
            record.update(diagnostics)
            record.update(n_treatment=int(samples["target"].eq("treatment").sum()),
                          n_control=int(samples["target"].eq("control").sum()),
                          n_blocks=int(samples[block].nunique()) if block else 0)
            if result is not None:
                for key in ("context", "stratum", "perturbation", "analysis_level", "guide"):
                    result[key] = record[key]
                results.append(result[RESULT_COLUMNS])
        record.update(status="skipped" if reason else "ok", reason=reason)
        summaries.append(record)

    if mode == "kolf_resample":
        for _, frame in pools.groupby("contrast_id", sort=True, observed=True):
            indices = pools.index.get_indexer(frame.index)
            target = frame.loc[frame["target"].ne(control), "target"].iloc[0]
            evaluate(frame, indices, frame["target"].eq(target), frame["target"].eq(control),
                     target, frame["analysis_level"].iloc[0])
    else:
        for context, frame in pools.groupby("context", sort=True, observed=True):
            indices = pools.index.get_indexer(frame.index)
            for target in sorted(set(frame["target"]) - {control}):
                evaluate(frame, indices, frame["target"].eq(target), frame["target"].eq(control),
                         target, "target")
                if settings["guide_level"]:
                    guides = sorted(frame.loc[frame["target"].eq(target), "guide"].unique())
                    for guide in guides:
                        evaluate(frame, indices, frame["guide"].eq(guide), frame["target"].eq(control),
                                 target, "guide", guide)
            groups = null_guide_sets(frame, control, stable_seed(settings["random_seed"], context, "null"),
                                     settings["ntc_null_groups"])
            for index, guides in enumerate(groups):
                ntc = frame["target"].eq(control)
                evaluate(frame, indices, ntc & frame["guide"].isin(guides),
                         ntc & ~frame["guide"].isin(guides), f"NTC-null-{index}", "ntc_null")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results = pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=RESULT_COLUMNS)
    results.to_csv(outdir / "differential_genes.tsv.gz", sep="\t", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(outdir / "contrast_qc.tsv.gz", sep="\t", index=False)
    guide_agreement(results, target_genes).to_csv(outdir / "guide_agreement.tsv.gz", sep="\t", index=False)
    write_json(outdir / "de_meta.json", {
        "dataset": dataset, "settings": settings, "mode": mode,
        "pydeseq2_version": importlib.metadata.version("pydeseq2"), "package_versions": versions(),
        "fdr_family": "genes within each target/guide/null contrast and context",
        "normalization": "DESeq2 poscounts size factors; no log-cell-count covariate",
        "paper_reproduction": False,
        "reproduction_scope": "KOLF mode reproduces cell-pool sampling; model/filter options remain explicitly configured.",
        "target_gene_mapping": symbol_col or "var_names",
        "selection": "all assigned cells; no Mixscape selection",
        "independent_filtering": True, "cooks_filter": True, "refit_cooks": True,
        "kolf_guide_level_note": "KOLF resampling currently tests pooled targets, not separate guides",
        "ntc_null_diagnostics": json.loads(
            summary.loc[summary["analysis_level"].eq("ntc_null")].to_json(orient="records")),
        "calibration_note": "Diagnostic control comparisons; few contrasts and reused controls do not establish FDR calibration.",
        "inference_note": ("Overlapping KOLF resamples can inflate significance; these are not biological replicates."
                           if mode == "kolf_resample" else "Replication follows the explicitly declared metadata mode."),
    })


if __name__ == "__main__":
    run(snakemake.input.counts, json.loads(snakemake.params.dataset),
        json.loads(snakemake.params.settings), snakemake.params.outdir)
