#!/usr/bin/env python3
"""
Diagnose AnnData file structure without modifying it
Uses backed mode to avoid loading full dataset into memory
"""
import scanpy as sc
import argparse
from pathlib import Path
import os

def diagnose_h5ad(input_h5ad, pert_col="gene"):
    """
    Diagnose AnnData file structure using backed mode
    """
    print("=" * 60)
    print("DIAGNOSING ANNDATA FILE")
    print("=" * 60)
    
    # File info
    print(f"\nFile: {input_h5ad}")
    if not os.path.exists(input_h5ad):
        print(f"ERROR: File does not exist!")
        return
    
    size_bytes = os.path.getsize(input_h5ad)
    size_gb = size_bytes / (1024**3)
    print(f"File size: {size_gb:.2f} GB")
    
    # Load in backed mode (doesn't load into memory)
    print("\nLoading file in backed mode...")
    adata = sc.read_h5ad(input_h5ad, backed='r')
    
    # Diagnose main matrix
    print("\n" + "=" * 60)
    print("MAIN MATRIX (.X)")
    print("=" * 60)
    print(f"Shape: {adata.shape}")
    print(f"Cells: {adata.n_obs:,}")
    print(f"Genes: {adata.n_vars:,}")
    print(f"Matrix type: {type(adata.X)}")
    
    # Check gene identifiers
    print(f"\nGene identifiers:")
    print(f"  var_names (first 5): {list(adata.var_names[:5])}")
    if 'gene_name' in adata.var.columns:
        print(f"  gene_name column (first 5): {list(adata.var['gene_name'][:5])}")
    
    print(f"\nVariable (.var) columns: {list(adata.var.columns)}")
    
    # Diagnose .raw layer
    print("\n" + "=" * 60)
    print("RAW LAYER (.raw)")
    print("=" * 60)
    if adata.raw is not None:
        print("Status: EXISTS")
        print(f"Shape: {adata.raw.shape}")
        print(f"Cells: {adata.raw.n_obs:,}")
        print(f"Genes: {adata.raw.n_vars:,}")
        
        # Check what's in raw
        print(f"\nRaw gene identifiers:")
        print(f"  var_names (first 5): {list(adata.raw.var_names[:5])}")
        if hasattr(adata.raw, 'var') and 'gene_name' in adata.raw.var.columns:
            print(f"  gene_name column (first 5): {list(adata.raw.var['gene_name'][:5])}")
        
        # Estimate raw size contribution
        raw_fraction = (adata.raw.n_vars / adata.n_vars) if adata.n_vars > 0 else 0
        estimated_raw_size = size_gb * 0.7 * raw_fraction  # Rough estimate
        print(f"\nEstimated .raw contribution: ~{estimated_raw_size:.2f} GB")
        
    else:
        print("Status: DOES NOT EXIST")
    
    # Check perturbation column
    print("\n" + "=" * 60)
    print(f"PERTURBATION COLUMN ('{pert_col}')")
    print("=" * 60)
    if pert_col in adata.obs.columns:
        unique_perts = adata.obs[pert_col].nunique()
        print(f"Unique perturbations: {unique_perts:,}")
        print(f"\nTop 10 most frequent perturbations:")
        print(adata.obs[pert_col].value_counts().head(10))
    else:
        print(f"WARNING: Column '{pert_col}' NOT FOUND")
        print(f"Available .obs columns: {list(adata.obs.columns)}")
    
    # Check for Mixscape results
    print("\n" + "=" * 60)
    print("MIXSCAPE ANNOTATIONS")
    print("=" * 60)
    mixscape_cols = [col for col in adata.obs.columns if 'mixscape' in col.lower()]
    if mixscape_cols:
        print(f"Found {len(mixscape_cols)} Mixscape columns:")
        for col in mixscape_cols:
            print(f"  - {col}")
            if adata.obs[col].dtype == 'object' or adata.obs[col].dtype.name == 'category':
                print(f"    Unique values: {adata.obs[col].nunique()}")
    else:
        print("No Mixscape columns found")
    
    # Summary and recommendation
    print("\n" + "=" * 60)
    print("SUMMARY & RECOMMENDATION")
    print("=" * 60)
    
    if adata.raw is not None:
        print(f"✗ File contains .raw layer with {adata.raw.n_vars:,} genes")
        print(f"✓ Main matrix has {adata.n_vars:,} genes")
        print(f"\nRECOMMENDATION:")
        print(f"  The .raw layer is taking up significant space.")
        print(f"  If you don't need it, consider removing it to:")
        print(f"  - Reduce file size by ~{size_gb * 0.5:.1f} GB")
        print(f"  - Speed up loading times")
        print(f"  - Simplify downstream processing")
    else:
        print(f"✓ File contains only main matrix with {adata.n_vars:,} genes")
        print(f"✓ No .raw layer present")
        print(f"\nFile is in optimal format for downstream analysis")
    
    # Close file
    adata.file.close()
    print("\n" + "=" * 60)
    print("Diagnosis complete!")
    print("=" * 60)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Diagnose AnnData file structure using backed mode (memory efficient)"
    )
    parser.add_argument("--input_h5ad", type=Path, required=True,
                        help="Input H5AD file to diagnose")
    parser.add_argument("--pert_col", type=str, default="gene",
                        help="Perturbation column name (default: gene)")
    
    args = parser.parse_args()
    
    diagnose_h5ad(args.input_h5ad, args.pert_col)