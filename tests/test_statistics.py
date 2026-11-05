"""Scientific contracts for count aggregation, replication and distribution tests."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from statistics_support import (
    aggregate, close_data, integer_counts, open_data, regular_membership,
    resample_membership, stable_seed,
)
from run_pseudobulk_de import collapse_biological_samples, contrast_data, guide_agreement
from run_distribution_tests import ad_test, energy_test, fit_basis


class Statistics(unittest.TestCase):
    def test_integer_counts_reject_fractional_and_negative_input(self):
        for values in ([[0.000001, 4]], [[-1, 4]], [[np.inf, 4]]):
            with self.assertRaises(ValueError):
                integer_counts(np.asarray(values))

    def test_streaming_counts_match_across_csr_csc_dense(self):
        values = np.arange(24).reshape(6, 4)
        obs = pd.DataFrame({"target": ["T", "T", "NT", "NT", "NT", "NT"],
                            "guide": ["g", "g", "n1", "n1", "n2", "n2"],
                            "sample": ["S"] * 6}, index=[f"c{i}" for i in range(6)])
        with tempfile.TemporaryDirectory() as tmp:
            for convert in (sparse.csr_matrix, sparse.csc_matrix, np.asarray):
                path = Path(tmp) / "input.h5ad"
                ad.AnnData(convert(values), obs=obs, var=pd.DataFrame(index=list("ABCD"))).write_h5ad(path)
                dataset = {"h5ad": str(path), "pert_col": "target", "control_label": "NT",
                           "statistics": {"sample_col": "sample", "guide_col": "guide"}}
                source, handle, meta = open_data(dataset)
                try:
                    membership, pools = regular_membership(meta, 6)
                    counts = aggregate(source, handle, membership, 2, 100).toarray()
                    for index, pool in pools.iterrows():
                        expected = values[obs["guide"].eq(pool["guide"])].sum(axis=0)
                        np.testing.assert_array_equal(counts[index], expected)
                    np.testing.assert_array_equal(np.asarray(membership.sum(axis=0)).ravel(), 1)
                finally:
                    close_data(source, handle)

    def test_biological_samples_collapse_guides_without_inventing_replicates(self):
        pools = pd.DataFrame({"context": ["ALL"] * 4, "target": ["T"] * 4,
                              "sample": ["S1", "S1", "S2", "S2"], "n_cells": [5] * 4})
        counts, samples = collapse_biological_samples(sparse.csr_matrix([[1], [2], [3], [4]]), pools)
        np.testing.assert_array_equal(counts.toarray(), [[3], [7]])
        self.assertEqual(samples["n_cells"].tolist(), [10, 10])

    def test_guide_calls_select_and_relabel_the_declared_cohort(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            obs = pd.DataFrame({"raw_target": ["wrong"] * 3, "guide": ["g", "unused", "ntc"],
                                "sample": ["S"] * 3}, index=["c0", "c1", "c2"])
            ad.AnnData(sparse.csr_matrix([[1], [100], [2]]), obs=obs,
                       var=pd.DataFrame(index=["GENE"])).write_h5ad(tmp / "counts.h5ad")
            pd.DataFrame({"cell_barcode": ["c2", "c0"], "gene_target": ["NT", "T"]}).to_csv(
                tmp / "calls.csv", index=False)
            dataset = {"h5ad": str(tmp / "counts.h5ad"), "guide_calls_csv": str(tmp / "calls.csv"),
                       "pert_col": "raw_target", "control_label": "NT",
                       "statistics": {"sample_col": "sample", "guide_col": "guide"}}
            source, handle, meta = open_data(dataset)
            try:
                self.assertEqual(meta["target"].tolist(), ["NT", "T"])
                membership, _ = regular_membership(meta, 3)
                self.assertEqual(aggregate(source, handle, membership, 2, 100).sum(), 3)
            finally:
                close_data(source, handle)

    def test_pool_size_filter_applies_after_biological_aggregation(self):
        pools = pd.DataFrame({
            "context": ["ALL"] * 12, "target": ["T"] * 6 + ["NT"] * 6,
            "sample": ["S1"] * 3 + ["S2"] * 3 + ["S1"] * 3 + ["S2"] * 3, "n_cells": [2] * 12})
        matrix, samples, reason = contrast_data(
            sparse.csr_matrix(np.ones((12, 1), dtype=int)), pools, pools["target"].eq("T"),
            pools["target"].eq("NT"), "biological_sample", True,
            {"min_cells_per_pool": 5, "min_replicates": 2})
        self.assertEqual(reason, "")
        self.assertEqual(matrix.shape, (4, 1))
        self.assertTrue(samples["n_cells"].eq(6).all())

    def test_one_batch_kolf_is_not_three_biological_samples(self):
        pools = pd.DataFrame({"context": ["ALL"] * 6, "target": ["T"] * 3 + ["NT"] * 3,
                              "sample": ["ALPHA"] * 6, "n_cells": [20] * 6})
        cfg = {"min_cells_per_pool": 5, "min_replicates": 3}
        matrix, _, reason = contrast_data(sparse.csr_matrix(np.ones((6, 10), dtype=int)),
                                          pools, pools["target"].eq("T"), pools["target"].eq("NT"),
                                          "biological_sample", True, cfg)
        self.assertIsNone(matrix)
        self.assertEqual(reason, "insufficient_matched_samples_or_pools")

    def test_resampling_overlaps_and_never_uses_target_cells_as_controls(self):
        meta = pd.DataFrame({"context": ["ALL"] * 40, "sample": ["S"] * 40, "batch": ["B"] * 40,
                             "target": ["T"] * 20 + ["NT"] * 20,
                             "guide": ["g"] * 20 + ["n"] * 20, "_row": np.arange(40)})
        membership, pools = resample_membership(meta, 40, "NT", 1, 3, 0.85)
        self.assertEqual(pools["n_cells"].tolist(), [17] * 6)
        self.assertTrue((np.asarray(membership[:, :20].sum(axis=0)).ravel() > 1).any())
        self.assertEqual(membership[pools["target"].eq("T").to_numpy(), 20:].nnz, 0)
        self.assertEqual(membership[pools["target"].eq("NT").to_numpy(), :20].nnz, 0)

    def test_seed_is_stable_across_processes(self):
        script = f"import sys; sys.path.insert(0, {str(ROOT / 'workflow/scripts')!r}); " \
                 "from statistics_support import stable_seed; print(stable_seed(3, 'T', 'B'))"
        result = subprocess.check_output([sys.executable, "-c", script], text=True)
        self.assertEqual(int(result), stable_seed(3, "T", "B"))

    def test_ad_feature_selection_does_not_require_pca(self):
        source = ad.AnnData(sparse.csr_matrix(np.arange(1, 31).reshape(6, 5)))
        settings = {"max_fit_cells": 6, "n_top_genes": 5, "pca_dims": 20,
                    "max_dense_elements": 1000, "normalize_target_sum": 10000}
        model, genes, _ = fit_basis(source, np.arange(6), settings, 7, need_pca=False)
        self.assertIsNone(model)
        # The middle column is constant after library normalization.
        self.assertEqual(set(genes), {0, 1, 3, 4})

    def test_guide_agreement_excludes_on_target_gene_using_explicit_mapping(self):
        rows = []
        for guide, values in (("a", [100, 1, 2, 3]), ("b", [100, 3, 2, 1])):
            for gene, value in zip(["ENSG1", "ENSG2", "ENSG3", "ENSG4"], values):
                rows.append(dict(context="ALL", perturbation="TP53", analysis_level="guide",
                                 guide=guide, gene=gene, log2FoldChange=value))
        results = pd.DataFrame(rows)
        agreement = guide_agreement(results, {"TP53": {"ENSG1"}})
        self.assertAlmostEqual(agreement["pearson_r"].item(), -1)
        unmapped = guide_agreement(results, {})
        self.assertTrue(pd.isna(unmapped["pearson_r"].item()))
        self.assertEqual(unmapped["reason"].item(), "unmapped_target_gene")

    def test_energy_matches_direct_formula_and_has_nonzero_permutation_pvalue(self):
        a, b = np.array([[0.], [1.], [2.]]), np.array([[4.], [5.], [6.]])
        direct = (2 * np.abs(a - b.T).mean() - np.abs(a - a.T).mean()
                  - np.abs(b - b.T).mean())
        distance, pvalue = energy_test(a, b, 99, 7, 100)
        self.assertAlmostEqual(distance, direct)
        self.assertGreaterEqual(pvalue, 1 / 100)
        self.assertEqual(energy_test(a, b, 99, 7, 100), (distance, pvalue))
        with self.assertRaisesRegex(ValueError, "budget"):
            energy_test(a, b, 99, 7, 1)

    def test_ad_detects_shape_change_at_equal_mean(self):
        a = np.tile([-2., 2.], 20)
        b = np.tile([-0.1, 0.1], 20)
        self.assertEqual(a.mean(), b.mean())
        _, pvalue = ad_test(a, b, 199, 7)
        self.assertLess(pvalue, 0.05)
        self.assertEqual(ad_test(np.zeros(20), np.zeros(20), 99, 7), (0., 1.))


if __name__ == "__main__":
    unittest.main()
