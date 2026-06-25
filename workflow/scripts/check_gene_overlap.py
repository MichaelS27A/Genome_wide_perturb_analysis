import scanpy as sc
import pandas as pd
import argparse
from pathlib import Path

def calculate_overlap(adata_path, genelist_path, adata_col, genelist_col):
    """
    Calculates the overlap between genes in an AnnData object and a gene list.
    """
    print("--- Calculating Gene Overlap ---")

    # 1. Load gene list of interest
    print(f"  - Loading your gene list from: {genelist_path}")
    your_genes_df = pd.read_csv(genelist_path, sep='\t')
    your_genes = set(your_genes_df[genelist_col].astype(str))
    print(f"  - Found {len(your_genes)} unique genes in your list.")

    # 2. Load AnnData in backed mode and get its genes
    print(f"  - Loading data genes from: {adata_path}")
    adata = sc.read_h5ad(adata_path, backed='r')
    data_genes = set(adata.var[adata_col].astype(str))
    adata.file.close()
    print(f"  - Found {len(data_genes)} unique genes in the data file.")

    # 3. Calculate overlap
    overlap = your_genes.intersection(data_genes)
    
    # 4. Print report
    print("\n--- Overlap Report ---")
    print(f"Genes in your list: {len(your_genes)}")
    print(f"Genes in the data file: {len(data_genes)}")
    print(f"Overlapping genes: {len(overlap)}")
    if len(your_genes) > 0:
        overlap_pct = (len(overlap) / len(your_genes)) * 100
        print(f"Overlap Percentage: {overlap_pct:.2f}% of your list is present in the data.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate gene list overlap with an AnnData file.")
    parser.add_argument("--input_h5ad", type=Path, required=True)
    parser.add_argument("--gene_list", type=Path, required=True)
    parser.add_argument("--adata_col", type=str, default="gene_name", help="Column in .var with gene symbols.")
    parser.add_argument("--genelist_col", type=str, default="GeneSymbol", help="Column in your gene list with symbols.")
    args = parser.parse_args()

    calculate_overlap(args.input_h5ad, args.gene_list, args.adata_col, args.genelist_col)