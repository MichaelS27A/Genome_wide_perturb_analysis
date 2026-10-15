"""Run the real classifier and DE workflow twice and verify numerical reproducibility."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "examples"))
from make_tiny_dataset import make_example


def run(config):
    subprocess.run([sys.executable, "-m", "snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
                    "--directory", str(ROOT), "--configfile", str(config),
                    "--cores", "1", "--rerun-incomplete"], check=True)


def verify(results):
    meta = json.loads((results / "run_provenance.json").read_text())
    assert list(meta["resolved_config"]["datasets"]) == ["TINY"], "Default atlases leaked into custom config"
    selected = pd.read_csv(results / "TINY/merged/selected_perturbed_cells_merged.tsv.gz", sep="\t")
    assert selected["cell_barcode"].is_unique and selected["perturbation"].nunique() == 4
    assert len(selected) > 40, "Synthetic strong perturbations should yield selected cells"
    clusters = pd.read_csv(results / "TINY/merged/perturbation_clusters.tsv", sep="\t")
    assert clusters.empty, "Different chunk PCA bases must not be clustered together"
    profiles = ad.read_h5ad(results / "TINY/postprocess/perturbation_pseudobulk.h5ad")
    assert profiles.shape == (4, 160)
    de = pd.read_csv(results / "TINY/global_de/perturbation_differential_genes_full.tsv.gz", sep="\t")
    self_de = de.loc[de["perturbation"] == de["gene"]]
    assert len(self_de) == 4 and (self_de["logfoldchange"] < -2).all()
    assert de[["perturbation", "gene"]].duplicated().sum() == 0
    counts = pd.read_csv(results / "TINY/global_de/group_cell_counts.tsv.gz", sep="\t")
    assert counts.loc[counts.iloc[:, 0] == "Non-Targeting", "n_cells"].item() == 60
    for path in (results / "TINY/chunk_runs").glob("*/run_meta.json"):
        chunk = json.loads(path.read_text())
        assert chunk["effective_parameters"]["pca_n_comps"] == 10
        assert chunk["package_versions"]["pertpy"] == "1.1.1"


def main(outdir):
    config = make_example(outdir)
    run(config)
    cfg = yaml.safe_load(config.read_text())
    first = Path(cfg["results_dir"])
    verify(first)
    cfg["results_dir"] = str(first.parent / "repeat_results")
    repeat = config.parent / "repeat.yaml"
    repeat.write_text(yaml.safe_dump(cfg, sort_keys=False))
    run(repeat)
    second = Path(cfg["results_dir"])
    verify(second)
    for relative in ("TINY/merged/selected_perturbed_cells_merged.tsv.gz",
                     "TINY/merged/perturbation_stats_merged.tsv.gz",
                     "TINY/postprocess/perturbation_umap_leiden.tsv.gz",
                     "TINY/global_de/perturbation_differential_genes_full.tsv.gz"):
        pd.testing.assert_frame_equal(pd.read_csv(first / relative, sep="\t"),
                                      pd.read_csv(second / relative, sep="\t"), rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(ad.read_h5ad(first / "TINY/postprocess/perturbation_pseudobulk.h5ad").X,
                               ad.read_h5ad(second / "TINY/postprocess/perturbation_pseudobulk.h5ad").X)
    print("PASS: real two-chunk workflow, input contracts and repeated numerical outputs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    main(args.outdir)
