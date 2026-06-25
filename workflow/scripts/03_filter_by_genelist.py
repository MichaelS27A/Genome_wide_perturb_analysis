import anndata as ad
import pandas as pd
import argparse
from pathlib import Path

def filter_h5ad(adata_path, gene_list_path, output_path, h5ad_symbol_col, genelist_symbol_col):
    """
    Filters an AnnData object to include only genes from a specified list.
    """
    print(f"--- Filtering AnnData file: {adata_path} ---")

    # 1. Load Gene List
    print(f"  - Loading gene list from: {gene_list_path}")
    try:
        gene_df = pd.read_csv(gene_list_path, sep='\t')
        if genelist_symbol_col not in gene_df.columns:
            raise ValueError(f"Column '{genelist_symbol_col}' not found in gene list file.")
        genes_to_keep = set(gene_df[genelist_symbol_col].unique())
        print(f"  - Found {len(genes_to_keep)} unique genes to keep.")
    except Exception as e:
        print(f"Error: Could not read or process gene list file. Reason: {e}")
        return

    # 2. Load AnnData and find matching genes
    print(f"  - Loading AnnData from: {adata_path}")
    try:
        adata = ad.read_h5ad(adata_path)
        # The primary index (.var_names) should contain the identifiers to match
        # This script assumes the .var_names hold the same type of identifier as your gene list.
        keep_mask = adata.var_names.isin(genes_to_keep)
        
        num_found = keep_mask.sum()
        if num_found == 0:
            print("Warning: No matching genes found between the list and the AnnData file. The output will be empty.")
        else:
            print(f"  - Found {num_found} matching genes in the AnnData file.")

    except Exception as e:
        print(f"Error: Could not read or process AnnData file. Reason: {e}")
        return

    # 3. Filter and Save
    filtered_adata = adata[:, keep_mask].copy()
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    filtered_adata.write_h5ad(output_path)
    print(f"--- Filtered AnnData with {filtered_adata.n_vars} genes saved to {output_path} ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter an .h5ad file based on a gene list.")
    parser.add_argument("--input_h5ad", type=Path, required=True)
    parser.add_argument("--gene_list", type=Path, required=True)
    parser.add_argument("--output_h5ad", type=Path, required=True)
    # These are kept for compatibility but the logic now assumes matching is on the index
    parser.add_argument("--h5ad_symbol_col", type=str, default="", help="Legacy. Not used.")
    parser.add_argument("--genelist_symbol_col", type=str, default="GeneSymbol", help="Column in gene list with symbols.")
    args = parser.parse_args()

    filter_h5ad(args.input_h5ad, args.gene_list, args.output_h5ad, args.h5ad_symbol_col, args.genelist_symbol_col)
