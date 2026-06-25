import scanpy as sc
import numpy as np
import argparse
from pathlib import Path

def qc_raw_data(adata_path, output_qc_path):
    """
    Performs a basic QC check on a raw AnnData object.
    """
    print(f"--- Performing QC on raw data file: {adata_path} ---")
    
    try:
        adata = sc.read_h5ad(adata_path, backed='r')
    except Exception as e:
        print(f"Error: Could not read file. Reason: {e}")
        return

    report_lines = []
    report_lines.append(f"QC Report for: {adata_path.name}")
    report_lines.append("="*40)
    report_lines.append(f"Shape: {adata.shape[0]} observations (cells) x {adata.shape[1]} variables (genes)")
    
    # Check if data looks like raw counts (simplified check)
    try:
        sample_data = adata.X[:10].toarray() if hasattr(adata.X, 'toarray') else adata.X[:10]
        max_val = sample_data.max()
        is_integer = np.allclose(sample_data, sample_data.astype(int), equal_nan=True)
        
        if max_val > 20 and is_integer:
            is_raw = "Looks like raw counts (integer data with reasonable range)."
        else:
            is_raw = "Data does not look like raw counts (float or low values)."
    except:
        is_raw = "Could not determine data type from sample."
        
    report_lines.append(f"Data type check: {is_raw}")
    
    report_lines.append("\n--- .obs columns ---")
    report_lines.append(str(adata.obs.head()))
    
    report_lines.append("\n--- .var columns ---")
    report_lines.append(str(adata.var.head()))
    
    adata.file.close()
    
    # Save report
    output_qc_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_qc_path, 'w') as f:
        f.write("\n".join(report_lines))
        
    print(f"--- QC Report saved to: {output_qc_path} ---")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform basic QC on a raw .h5ad file.")
    parser.add_argument("--input_h5ad", type=Path, required=True)
    parser.add_argument("--output_qc", type=Path, required=True)
    args = parser.parse_args()
    qc_raw_data(args.input_h5ad, args.output_qc)
