"""Shared H5AD access and strict table contracts for the chunk workflow."""

import os
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import hashlib
import importlib.metadata
import platform
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, sparse_dataset

def _open_minimal_backed_adata(h5ad_path: Path) -> tuple[ad.AnnData, h5py.File | None]:
    """Open H5AD without eagerly materializing heavy optional groups (e.g. layers).

    Falls back to anndata backed reader for files whose X metadata does not expose
    the sparse encoding attributes expected by sparse_dataset().
    """
    h5 = h5py.File(h5ad_path, "r")
    x_obj = h5.get("X")
    has_sparse_attrs = (
        isinstance(x_obj, h5py.Group)
        and ("encoding-type" in x_obj.attrs or "h5sparse_format" in x_obj.attrs)
    )
    if not has_sparse_attrs and not isinstance(x_obj, h5py.Dataset):
        h5.close()
        return ad.read_h5ad(h5ad_path, backed="r"), None
    try:
        X = x_obj if isinstance(x_obj, h5py.Dataset) else sparse_dataset(x_obj)
        obs = read_elem(h5["obs"])
        var = read_elem(h5["var"])
        if not obs.index.is_unique or not var.index.is_unique:
            raise ValueError("H5AD cell and gene identifiers must be unique")
        adata = ad.AnnData(X=X, obs=obs, var=var)
        return adata, h5
    except Exception as e:
        h5.close()
        msg = str(e)
        if "encoding-type" not in msg and "h5sparse_format" not in msg:
            raise
        # Compatibility path for variant/legacy H5AD matrix encoding.
        adata = ad.read_h5ad(h5ad_path, backed="r")
        return adata, None


def _close_backing_resources(adata: ad.AnnData, h5_handle: h5py.File | None) -> None:
    if h5_handle is not None:
        try:
            h5_handle.close()
        except Exception:
            pass
        return
    file_handle = getattr(adata, "file", None)
    if file_handle is not None:
        try:
            file_handle.close()
        except Exception:
            pass


def _read_csc_nnz_per_gene(h5: h5py.File | None) -> tuple[np.ndarray | None, int]:
    """Read per-gene nnz for CSC X if available; returns (nnz_per_gene, bytes_per_nnz)."""
    if h5 is None:
        return None, 12
    try:
        xg = h5["X"]
        if "indptr" not in xg:
            return None, 12
        indptr = np.asarray(xg["indptr"][...], dtype=np.int64)
        if indptr.ndim != 1 or indptr.size < 2:
            return None, 12
        nnz = np.diff(indptr).astype(np.int64, copy=False)
        bytes_per_nnz = int(xg["data"].dtype.itemsize + xg["indices"].dtype.itemsize)
        if bytes_per_nnz <= 0:
            bytes_per_nnz = 12
        return nnz, bytes_per_nnz
    except Exception:
        return None, 12


def _select_csc_gene_subset(
    var: pd.DataFrame,
    nnz_per_gene: np.ndarray | None,
    max_genes: int,
    max_total_nnz: int,
) -> tuple[np.ndarray, dict[str, int]]:
    if max_genes <= 0 or max_genes >= int(var.shape[0]):
        full = np.arange(int(var.shape[0]), dtype=np.int64)
        return full, {
            "selected_genes": int(full.size),
            "selected_total_nnz": int(-1),
            "per_gene_nnz_cap": int(-1),
            "selector_used_budget": int(0),
        }

    n_vars = int(var.shape[0])
    if nnz_per_gene is None or int(nnz_per_gene.shape[0]) != n_vars:
        nnz_per_gene = np.zeros(n_vars, dtype=np.int64)
    nnz_per_gene = nnz_per_gene.astype(np.int64, copy=False)

    score_cols = (
        "highly_variable_rank",
        "dispersions_norm",
        "variance",
        "mean_counts",
        "total_counts",
        "n_cells_by_counts",
        "n_cells",
    )
    order = np.arange(n_vars, dtype=np.int64)
    for col in score_cols:
        if col not in var.columns:
            continue
        scores = pd.to_numeric(var[col], errors="coerce").to_numpy(dtype=np.float64, copy=False)
        valid = np.isfinite(scores)
        if not bool(valid.any()):
            continue
        if col == "highly_variable_rank":
            safe_scores = np.where(valid, scores, np.inf)
            order = np.argsort(safe_scores, kind="stable").astype(np.int64, copy=False)
        else:
            safe_scores = np.where(valid, scores, -np.inf)
            order = np.argsort(safe_scores, kind="stable")[::-1].astype(np.int64, copy=False)
        break

    nonzero_mask = nnz_per_gene > 0
    selected: list[int] = []
    selected_set: set[int] = set()
    selected_total_nnz = 0

    if max_total_nnz > 0:
        per_gene_cap = max(1, int(max_total_nnz) // max(1, int(max_genes)))
        min_candidates = max(100, int(max_genes) // 3)
        candidate_order = np.array([], dtype=np.int64)
        for fac in (1, 2, 4, 8, 16, 32, 64, 128):
            cap = int(per_gene_cap) * int(fac)
            cand_mask = nonzero_mask & (nnz_per_gene <= cap)
            candidate_order = order[cand_mask[order]]
            if int(candidate_order.size) >= min_candidates or fac == 128:
                per_gene_cap = cap
                break
    else:
        per_gene_cap = -1
        candidate_order = order[nonzero_mask[order]]

    for g in candidate_order:
        gi = int(g)
        gn = int(nnz_per_gene[gi])
        if max_total_nnz > 0 and selected_total_nnz + gn > int(max_total_nnz):
            continue
        selected.append(gi)
        selected_set.add(gi)
        selected_total_nnz += gn
        if len(selected) >= int(max_genes):
            break

    if len(selected) < int(max_genes):
        remainder = order[
            nonzero_mask[order]
            & np.array([int(x) not in selected_set for x in order], dtype=bool)
        ]
        for g in remainder:
            gi = int(g)
            gn = int(nnz_per_gene[gi])
            if max_total_nnz > 0 and selected_total_nnz + gn > int(max_total_nnz):
                continue
            selected.append(gi)
            selected_set.add(gi)
            selected_total_nnz += gn
            if len(selected) >= int(max_genes):
                break

    if not selected:
        raise ValueError("CSC feature budget cannot select any nonzero genes; increase the explicit budget")

    out = np.sort(np.asarray(selected, dtype=np.int64))
    return out, {
        "selected_genes": int(out.size),
        "selected_total_nnz": int(selected_total_nnz),
        "per_gene_nnz_cap": int(per_gene_cap),
        "selector_used_budget": int(max_total_nnz > 0),
    }


def require_unique(frame, keys, context):
    missing = set(keys) - set(frame.columns)
    if missing:
        raise ValueError(f"{context}: missing columns {sorted(missing)}")
    if frame[keys].isna().any().any():
        raise ValueError(f"{context}: missing identifiers or labels")
    if frame.duplicated(keys).any():
        raise ValueError(f"{context}: duplicate identifiers {keys}")


def require_barcodes(obs_names, barcodes):
    if not obs_names.is_unique:
        raise ValueError("H5AD observation names must be unique")
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("Chunk cell barcodes must be unique")
    missing = pd.Index(barcodes).difference(obs_names)
    if len(missing):
        raise ValueError(f"{len(missing)} chunk barcodes missing from H5AD: {missing[:5].tolist()}")


def read_tables(paths):
    # A header-only table is valid; a missing or corrupt chunk is a failure.
    return [pd.read_csv(p, sep="\t", compression="infer") for p in paths]


def axis_hash(values):
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def runtime_metadata(adata):
    return {
        "python": platform.python_version(),
        "package_versions": {name: importlib.metadata.version(name) for name in
            ("pertpy", "scanpy", "anndata", "numpy", "pandas", "scipy", "scikit-learn", "h5py")},
        "cell_order_sha256": axis_hash(adata.obs_names),
        "gene_order_sha256": axis_hash(adata.var_names),
    }


def require_raw_counts(matrix):
    from scipy import sparse
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("H5AD X must contain finite, nonnegative raw integer counts")
