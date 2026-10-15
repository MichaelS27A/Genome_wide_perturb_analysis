"""Generate a deterministic raw-count screen; use only to verify software."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import yaml


def make_example(outdir, label_source="csv", storage="csr"):
    outdir = Path(outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(27)
    labels = np.repeat(["Non-Targeting", "TARGET_A", "TARGET_B", "TARGET_C", "TARGET_D"], 80)
    means = np.tile(rng.uniform(1, 5, 160), (len(labels), 1))
    # Independent target-specific expression programs with clear knockdown.
    for i, target in enumerate(("TARGET_A", "TARGET_B", "TARGET_C", "TARGET_D")):
        mask = labels == target
        means[mask, i] = 0.05
        means[mask, 10 + i * 25:35 + i * 25] *= 4
    counts = rng.poisson(means).astype(np.float32)
    barcodes = [f"cell_{i:04d}" for i in range(len(labels))]
    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    if label_source == "obs":
        obs["gene_target"] = labels
    var = pd.DataFrame(index=["TARGET_A", "TARGET_B", "TARGET_C", "TARGET_D"] + [f"GENE_{i}" for i in range(156)])
    matrix = {"csr": sparse.csr_matrix, "csc": sparse.csc_matrix, "dense": np.asarray}[storage](counts)
    ad.AnnData(matrix, obs=obs, var=var).write_h5ad(outdir / "counts.h5ad")
    calls = pd.DataFrame({"cell_barcode": barcodes, "gene_target": labels})
    calls.to_csv(outdir / "calls.csv", index=False)
    repo = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((repo / "config/config.example.yaml").read_text())
    config["datasets"] = {"TINY": {
        "h5ad": str(outdir / "counts.h5ad"), "pert_col": "gene_target",
        "control_label": "Non-Targeting", "data_id": "synthetic-software-test-v1",
    }}
    if label_source == "csv":
        config["datasets"]["TINY"]["guide_calls_csv"] = str(outdir / "calls.csv")
    config["results_dir"] = str(outdir / "results")
    config["chunking"].update(perturbations_per_chunk=2, min_cells_per_perturbation=5,
                              max_controls_per_chunk=60, csv_read_chunk_size=57)
    config["mixscape"].update(pca_n_comps=10, pca_dims=10, signature_n_dims=10,
                              signature_n_neighbors=10, min_de_genes=5, mem_mb=4000)
    config["postprocess"].update(min_selected_cells=5, n_top_de_genes=20)
    config["global_de"] = {"n_top_de_genes": 20}
    config["methods"]["global_de"]["enabled"] = True
    (outdir / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    return outdir / "config.yaml"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("examples/tiny"))
    parser.add_argument("--label-source", choices=("csv", "obs"), default="csv")
    parser.add_argument("--storage", choices=("csr", "csc", "dense"), default="csr")
    args = parser.parse_args()
    print(make_example(args.outdir, args.label_source, args.storage))
