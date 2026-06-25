import pandas as pd
import argparse
from pathlib import Path

def prepare_lists(annotation_file, output_dir, gene_col, function_col):
    """
    Reads a comprehensive gene annotation file and creates separate lists
    of TFs and target genes (transporters/enzymes).
    """
    print(f"--- Preparing Gene Lists from Annotation File: {annotation_file} ---")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        df = pd.read_csv(annotation_file, sep='\t')
        print(f"  - Loaded annotation file with {len(df)} rows.")
    except Exception as e:
        print(f"Error: Could not read annotation file. Reason: {e}")
        return

    if gene_col not in df.columns or function_col not in df.columns:
        print(f"Error: Annotation file must contain '{gene_col}' and '{function_col}' columns.")
        return

    # --- Extract TFs ---
    # Find rows where the function column contains 'TF' or 'transcription factor'
    tf_mask = df[function_col].str.contains("TF|transcription factor", case=False, na=False)
    tfs = df[tf_mask][gene_col].unique()
    tf_df = pd.DataFrame(tfs)
    tf_output_path = output_dir / "tfs.txt"
    tf_df.to_csv(tf_output_path, header=False, index=False)
    print(f"  - Found {len(tfs)} TFs. Saved to: {tf_output_path}")

    # --- Extract Targets ---
    # Find rows where the function column contains 'Transporter' or 'Enzyme'
    target_mask = df[function_col].str.contains("Transporter|Enzyme", case=False, na=False)
    targets = df[target_mask][gene_col].unique()
    target_df = pd.DataFrame(targets)
    target_output_path = output_dir / "targets.txt"
    target_df.to_csv(target_output_path, header=False, index=False)
    print(f"  - Found {len(targets)} targets (Transporters/Enzymes). Saved to: {target_output_path}")

    print("--- Gene List Preparation Complete ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare TF and target gene lists from an annotation file.")
    parser.add_argument("--annotation_file", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    parser.add_argument("--gene_col", type=str, default="GeneSymbol", help="Column name for gene symbols.")
    parser.add_argument("--function_col", type=str, default="GeneType", help="Column name for gene functions.")
    args = parser.parse_args()

    prepare_lists(args.annotation_file, args.output_dir, args.gene_col, args.function_col)
