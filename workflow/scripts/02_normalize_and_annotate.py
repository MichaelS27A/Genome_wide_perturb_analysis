#!/usr/bin/env python

import pertpy as pt
import scanpy as sc
import argparse
from pathlib import Path

def run_mixscape_analysis(input_h5ad, output_h5ad, pert_col, control_val, use_hvg_for_pca=False):
    """
    Runs the full Mixscape workflow including preprocessing and analysis.
    """
    print("--- Starting Mixscape Preprocessing and Analysis ---")

    # 1. Load Data
    print(f"  - Loading AnnData from: {input_h5ad}")
    try:
        adata = sc.read_h5ad(input_h5ad)
    except Exception as e:
        print(f"Error: Could not read input file. Reason: {e}")
        return

    # 2. Preprocessing (as described in the pertpy tutorial)
    print("  - Normalizing total counts and log-transforming...")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    print("  - Finding highly variable genes...")
    sc.pp.highly_variable_genes(adata, subset=False)

    print("  - Running PCA...")
    sc.pp.pca(adata, use_highly_variable=use_hvg_for_pca)

    # 3. Run Mixscape Workflow
    print("  - Initializing Mixscape...")
    ms = pt.tl.Mixscape()

    print("  - Calculating perturbation signatures...")
    # Using 'replicate' as split_by is common, but we make it optional.
    # The script will work even if the column doesn't exist.
    split_by_col = 'replicate' if 'replicate' in adata.obs.columns else None
    ms.perturbation_signature(adata, pert_key=pert_col, control=control_val, split_by=split_by_col)

    print("  - Running Mixscape classification...")
    ms.mixscape(adata=adata, control=control_val, labels=pert_col, layer="X_pert")

    print("--- Mixscape Classification Summary ---")
    print(adata.obs["mixscape_class_global"].value_counts())

    # 4. Save Annotated Data
    print(f"  - Saving annotated data to: {output_h5ad}")
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)

    print("--- Workflow Complete ---")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a preprocessing workflow including Mixscape.")
    parser.add_argument("--input_h5ad", type=Path, required=True)
    parser.add_argument("--output_h5ad", type=Path, required=True)
    parser.add_argument("--pert_col", type=str, required=True)
    parser.add_argument("--control_val", type=str, required=True)
    parser.add_argument(
        "--use_hvg_for_pca",
        action="store_true",
        help="If set, PCA is computed on highly variable genes only. Default uses full transcriptome.",
    )
    args = parser.parse_args()

    run_mixscape_analysis(
        args.input_h5ad,
        args.output_h5ad,
        args.pert_col,
        args.control_val,
        use_hvg_for_pca=args.use_hvg_for_pca,
    )
