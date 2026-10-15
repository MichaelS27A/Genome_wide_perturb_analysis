"""Regression checks for scientific input and merge integrity."""

import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml
from jsonschema import Draft202012Validator

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
from workflow_support import require_raw_counts, read_tables, _select_csc_gene_subset


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Contracts(unittest.TestCase):
    def test_config_rejects_typos_and_invalid_bounds(self):
        schema = yaml.safe_load((ROOT / "workflow/schemas/config.schema.yaml").read_text())
        config = yaml.safe_load((ROOT / "config/config.example.yaml").read_text())
        validator = Draft202012Validator(schema)
        validator.validate(config)
        config["mixscape"]["signature_n_neighbors"] = 0
        self.assertTrue(list(validator.iter_errors(config)))
        config["mixscape"]["signature_n_neighbors"] = 30
        config["mixscape"]["signature_n_neighbours"] = 30
        self.assertTrue(list(validator.iter_errors(config)))

    def test_matrix_storage_preserves_cells_and_raw_counts(self):
        runner = load("05_run_mixscape_chunk.py")
        values = np.arange(30).reshape(6, 5).astype(float)
        obs = pd.DataFrame(index=[f"c{i}" for i in range(6)])
        var = pd.DataFrame(index=[f"g{i}" for i in range(5)])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            calls = tmp / "cells.tsv"
            pd.DataFrame({"cell_barcode": ["c4", "c1"], "gene_target": ["T", "NT"],
                          "is_control": [0, 1]}).to_csv(calls, sep="\t", index=False)
            for convert in (sparse.csr_matrix, sparse.csc_matrix, np.asarray):
                path = tmp / "counts.h5ad"
                source = ad.AnnData(convert(values), obs=obs, var=var)
                # Optional layers must not silently replace the declared count X.
                source.layers["counts"] = convert(np.zeros_like(values))
                source.write_h5ad(path)
                subset = runner.load_subset(path, calls, "target", 0, 0)
                actual = subset.X.toarray() if sparse.issparse(subset.X) else subset.X
                np.testing.assert_array_equal(actual, values[[1, 4]])
                self.assertEqual(subset.obs["target"].tolist(), ["NT", "T"])
                pd.DataFrame({"cell_barcode": ["absent"], "gene_target": ["T"],
                              "is_control": [0]}).to_csv(calls, sep="\t", index=False)
                with self.assertRaisesRegex(ValueError, "missing from H5AD"):
                    runner.load_subset(path, calls, "target", 0, 0)
                pd.DataFrame({"cell_barcode": ["c4", "c1"], "gene_target": ["T", "NT"],
                              "is_control": [0, 1]}).to_csv(calls, sep="\t", index=False)

    def test_missing_chunk_is_not_an_empty_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            empty = tmp / "empty.tsv"
            pd.DataFrame(columns=["perturbation", "gene"]).to_csv(empty, sep="\t", index=False)
            self.assertEqual(read_tables([empty])[0].columns.tolist(), ["perturbation", "gene"])
            with self.assertRaises(FileNotFoundError):
                read_tables([empty, tmp / "missing.tsv"])

    def test_configured_classification_and_missing_output(self):
        runner = load("05_run_mixscape_chunk.py")
        obs = pd.DataFrame({"mixscape_class_global": ["KD", "NP", "NT"]})
        self.assertEqual(runner.detect_perturbed_mask(obs, "KD").tolist(), [True, False, False])
        with self.assertRaisesRegex(RuntimeError, "did not produce"):
            runner.detect_perturbed_mask(pd.DataFrame(index=[0]))

    def test_impossible_feature_budget_fails_instead_of_substituting_genes(self):
        with self.assertRaisesRegex(ValueError, "cannot select"):
            _select_csc_gene_subset(pd.DataFrame(index=["G1", "G2", "G3"]),
                                    np.array([100, 100, 100]), 2, 1)

    def test_duplicate_perturbations_and_mismatched_features_fail_merge(self):
        merger = load("20_merge_postprocess_chunk_results.py")
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp) / "a.h5ad", Path(tmp) / "b.h5ad"]
            for path in paths:
                ad.AnnData(np.ones((1, 2)), obs=pd.DataFrame(index=["T"]),
                           var=pd.DataFrame(index=["G1", "G2"])).write_h5ad(path)
            with self.assertRaisesRegex(ValueError, "Duplicate perturbations"):
                merger._merge_pseudobulk(paths)
            ad.AnnData(np.ones((1, 2)), obs=pd.DataFrame(index=["U"]),
                       var=pd.DataFrame(index=["G2", "G1"])).write_h5ad(paths[1])
            with self.assertRaisesRegex(ValueError, "different gene"):
                merger._merge_pseudobulk(paths)

    def test_small_dataset_has_no_fabricated_embedding(self):
        merger = load("20_merge_postprocess_chunk_results.py")
        pb = ad.AnnData(np.ones((2, 3)), obs=pd.DataFrame(index=["T", "U"]))
        result = merger._embed_merged_pseudobulk(pb)
        self.assertTrue(result[["umap_1", "umap_2"]].isna().all().all())
        self.assertTrue(result["leiden"].isna().all())

    def test_normalized_expression_cannot_be_used_as_counts(self):
        require_raw_counts(sparse.csr_matrix([[0, 4], [1, 0]]))
        for values in ([[0.5, 1]], [[-1, 2]], [[np.nan, 1]]):
            with self.assertRaisesRegex(ValueError, "raw integer counts"):
                require_raw_counts(np.asarray(values))

    def test_duplicate_guide_calls_across_streaming_blocks_fail(self):
        builder = load("04_build_mixscape_chunk_manifests.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            calls = tmp / "calls.csv"
            pd.DataFrame({"cell_barcode": ["c1", "c2", "c1"],
                          "gene_target": ["NT", "T", "T"]}).to_csv(calls, index=False)
            with self.assertRaisesRegex(ValueError, "duplicate cell barcodes"):
                builder.build_from_csv(calls, tmp, "NT", 2, 0, 1, 0, 1, 0)


if __name__ == "__main__":
    unittest.main()
