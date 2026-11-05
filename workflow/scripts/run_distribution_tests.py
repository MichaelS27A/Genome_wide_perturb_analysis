"""Seeded energy-distance and explicitly permutation-calibrated Anderson–Darling tests."""

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).parent))
from statistics_support import (
    bh, close_data, integer_counts, null_guide_sets, open_data, stable_seed,
    versions, write_json,
)

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import PermutationMethod, anderson_ksamp
from sklearn.decomposition import PCA


ENERGY_COLUMNS = ["context", "stratum", "perturbation", "analysis_level",
                  "n_treatment", "n_control", "energy_distance", "pvalue", "padj"]
AD_COLUMNS = ["context", "stratum", "perturbation", "analysis_level", "gene",
              "stat", "log2_mean_normalized_ratio", "pvalue", "padj"]


def energy_test(first, second, permutations, seed, max_elements):
    """Biased/V-statistic estimator, including within-population diagonal zeros."""
    joined = np.vstack([first, second])
    n, m = len(first), len(second)
    if min(n, m) < 2:
        raise ValueError("Energy testing needs at least two cells per population")
    if len(joined)**2 > max_elements:
        raise ValueError("Energy-distance matrix exceeds the declared element budget")
    distances = cdist(joined, joined, metric="euclidean")
    row_sums, total = distances.sum(axis=1), distances.sum()

    def statistic(weights):
        within_first = ((weights @ distances) * weights).sum(axis=1)
        first_total = weights @ row_sums
        within_second = total - 2 * first_total + within_first
        cross = first_total - within_first
        return 2 * cross / (n * m) - within_first / n**2 - within_second / m**2

    observed = statistic(np.r_[np.ones(n), np.zeros(m)][None, :])[0]
    rng = np.random.default_rng(seed)
    exceedances = 0
    for start in range(0, permutations, 64):
        weights = np.zeros((min(64, permutations - start), n + m))
        for row in weights:
            row[rng.choice(n + m, n, replace=False)] = 1
        exceedances += int((statistic(weights) >= observed - 1e-12).sum())
    return float(observed), (exceedances + 1) / (permutations + 1)


def ad_test(first, second, permutations, seed):
    if np.ptp(np.concatenate([first, second])) == 0:
        return 0.0, 1.0
    result = anderson_ksamp(
        [first, second],
        method=PermutationMethod(n_resamples=permutations, rng=np.random.default_rng(seed)),
    )
    return float(result.statistic), float(result.pvalue)


def normalized_rows(adata, rows, settings):
    rows = np.sort(np.asarray(rows, dtype=int))
    if not sparse.issparse(adata.X) and not hasattr(adata.X, "format"):
        if len(rows) * adata.n_vars > settings["max_dense_elements"]:
            raise ValueError("Raw dense subset exceeds the declared element budget")
    values = sparse.csr_matrix(integer_counts(adata.X[rows, :]), dtype=float)
    totals = np.asarray(values.sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise ValueError("Zero-library cells must be removed in the declared input QC")
    values = sparse.diags(settings["normalize_target_sum"] / totals) @ values
    values.data = np.log1p(values.data)
    return values.tocsr(), rows


def fit_basis(adata, rows, settings, seed, need_pca=True):
    rng = np.random.default_rng(seed)
    selected = rng.choice(rows, size=min(len(rows), settings["max_fit_cells"]), replace=False)
    values, selected = normalized_rows(adata, selected, settings)
    means = np.asarray(values.mean(axis=0)).ravel()
    variances = np.asarray(values.power(2).mean(axis=0)).ravel() - means**2
    candidates = np.flatnonzero(variances > 1e-12)
    order = candidates[np.argsort(-variances[candidates], kind="stable")]
    genes = np.sort(order[:settings["n_top_genes"]]) if settings["n_top_genes"] else np.sort(order)
    if not len(genes):
        raise ValueError("No variable genes in the declared fitting cohort")
    if need_pca and min(len(selected), len(genes)) <= settings["pca_dims"]:
        raise ValueError("Requested PCA dimensions exceed the available variable features/cells")
    if len(selected) * len(genes) > settings["max_dense_elements"]:
        raise ValueError("PCA fitting exceeds the declared dense element budget")
    dense = values[:, genes].toarray()
    # sklearn's random_state is 32-bit; the full sampling seed is recorded separately.
    model = (PCA(n_components=settings["pca_dims"], svd_solver="randomized",
                 random_state=int(seed % 2**32)).fit(dense) if need_pca else None)
    return model, genes, selected


def run(dataset, settings, energy_enabled, ad_enabled, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if settings["max_control_cells"] < settings["min_cells_per_group"]:
        raise ValueError("max_control_cells must be at least min_cells_per_group")
    adata, handle, meta = open_data(dataset)
    energy_rows, ad_rows, qc, basis_arrays, basis_meta = [], [], [], {}, []
    try:
        for stratum, frame in meta.groupby(["context", "sample", "batch"], sort=True):
            identity = json.dumps(list(stratum))
            seed = stable_seed(settings["random_seed"], identity, "basis")
            model, genes, fit_rows = fit_basis(adata, frame["_row"].to_numpy(), settings, seed, energy_enabled)
            key = f"S{len(basis_meta)}"
            basis_arrays.update({f"{key}_rows": fit_rows, f"{key}_genes": genes,
                                 f"{key}_components": model.components_ if model else np.empty((0, len(genes))),
                                 f"{key}_mean": model.mean_ if model else np.empty(0)})
            basis_meta.append(dict(key=key, stratum=identity, seed=str(seed),
                                   genes=adata.var_names[genes].tolist()))
            control = dataset["control_label"]
            ntc = frame["target"].eq(control)
            control_rows = frame.loc[ntc, "_row"].to_numpy()
            if len(control_rows) < settings["min_cells_per_group"]:
                raise ValueError(f"Insufficient matched controls for distribution stratum {identity}")
            tasks = [(target, cells["_row"].to_numpy(), control_rows, "target")
                     for target, cells in frame.loc[~ntc].groupby("target", sort=True)]
            for index, guides in enumerate(null_guide_sets(
                    frame, control, stable_seed(settings["random_seed"], identity, "null"),
                    settings["ntc_null_groups"])):
                held_out = ntc & frame["guide"].isin(guides)
                remaining = ntc & ~frame["guide"].isin(guides)
                tasks.append((f"NTC-null-{index}", frame.loc[held_out, "_row"].to_numpy(),
                              frame.loc[remaining, "_row"].to_numpy(), "ntc_null"))
            for target, rows, references, level in tasks:
                info = dict(context=stratum[0], stratum=identity, perturbation=target, analysis_level=level)
                if min(len(rows), len(references)) < settings["min_cells_per_group"]:
                    qc.append({**info, "status": "skipped", "reason": "insufficient_cells"})
                    continue
                rng = np.random.default_rng(stable_seed(settings["random_seed"], identity, "controls", level))
                references = rng.choice(references, min(len(references), settings["max_control_cells"]), replace=False)
                if max(len(rows), len(references)) * len(genes) > settings["max_dense_elements"]:
                    raise ValueError("Distribution contrast exceeds the declared dense element budget")
                first, _ = normalized_rows(adata, rows, settings)
                second, _ = normalized_rows(adata, references, settings)
                first, second = first[:, genes].toarray(), second[:, genes].toarray()
                test_seed = stable_seed(settings["random_seed"], identity, target, level)
                if energy_enabled:
                    distance, pvalue = energy_test(
                        model.transform(first), model.transform(second), settings["n_permutations"],
                        test_seed, settings["max_distance_elements"])
                    energy_rows.append({**info, "n_treatment": len(first), "n_control": len(second),
                                        "energy_distance": distance, "pvalue": pvalue})
                if ad_enabled:
                    effect = np.log2((np.expm1(first).mean(axis=0) + 0.5)
                                     / (np.expm1(second).mean(axis=0) + 0.5))
                    for index, gene in enumerate(genes):
                        statistic, pvalue = ad_test(first[:, index], second[:, index],
                                                   settings["n_permutations"],
                                                   stable_seed(test_seed, adata.var_names[gene]))
                        ad_rows.append({**info, "gene": str(adata.var_names[gene]), "stat": statistic,
                                        "log2_mean_normalized_ratio": float(effect[index]), "pvalue": pvalue})
                qc.append({**info, "status": "ok", "reason": "",
                           "n_treatment": len(rows), "n_control": len(references)})
    finally:
        close_data(adata, handle)
    energy = pd.DataFrame(energy_rows, columns=ENERGY_COLUMNS)
    genes = pd.DataFrame(ad_rows, columns=AD_COLUMNS)
    for frame, keys in ((energy, ["context", "analysis_level"]),
                        (genes, ["context", "stratum", "perturbation", "analysis_level"])):
        for _, indices in frame.groupby(keys, sort=False).groups.items():
            frame.loc[indices, "padj"] = bh(frame.loc[indices, "pvalue"])
    energy.to_csv(outdir / "energy_distance.tsv.gz", sep="\t", index=False)
    genes.to_csv(outdir / "anderson_darling_permutation.tsv.gz", sep="\t", index=False)
    pd.DataFrame(qc).to_csv(outdir / "contrast_qc.tsv.gz", sep="\t", index=False)
    np.savez_compressed(outdir / "basis.npz", **basis_arrays)
    write_json(outdir / "distribution_meta.json", {
        "dataset": dataset, "settings": settings, "package_versions": versions(), "basis": basis_meta,
        "energy_enabled": energy_enabled, "anderson_darling_permutation_enabled": ad_enabled,
        "preprocessing": "log1p library-normalized counts; variance-selected features on a fixed seeded fit sample",
        "paper_reproduction": False, "experimental_unit": "cells within declared context/sample/batch strata",
        "energy_estimator": "biased V statistic; Euclidean PCA distance; diagonal zeros included",
        "energy_fdr_family": "perturbation/stratum tests within each context; target and null families separate",
        "ad_fdr_family": "selected genes within each perturbation/stratum; permutation calibration",
        "minimum_monte_carlo_pvalue": 1 / (settings["n_permutations"] + 1),
        "effect_pseudocount": 0.5,
    })


if __name__ == "__main__":
    run(json.loads(snakemake.params.dataset), json.loads(snakemake.params.settings),
        snakemake.params.energy_enabled, snakemake.params.ad_enabled, snakemake.params.outdir)
