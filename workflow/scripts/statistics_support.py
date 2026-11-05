"""Count, metadata and sampling contracts shared by optional statistical analyses."""

import hashlib
import importlib.metadata
import json
import math
from pathlib import Path

from workflow_support import (
    _close_backing_resources, _open_minimal_backed_adata, require_barcodes,
    require_unique,
)

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control


def stable_seed(seed, *labels):
    text = json.dumps([int(seed), *map(str, labels)], ensure_ascii=False)
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def versions():
    return {name: importlib.metadata.version(name) for name in
            ("numpy", "pandas", "scipy", "anndata", "scikit-learn")}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def bh(pvalues):
    values = np.asarray(pvalues, dtype=float)
    valid = np.isfinite(values)
    result = np.full(values.shape, np.nan)
    if valid.any():
        if ((values[valid] < 0) | (values[valid] > 1)).any():
            raise ValueError("P-values must be in [0, 1]")
        result[valid] = false_discovery_control(values[valid], method="bh")
    return result


def integer_counts(matrix):
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if (not np.isfinite(values).all() or (values < 0).any()
            or (values != np.rint(values)).any() or (values > 2**53).any()):
        raise ValueError("Statistics require non-negative, finite raw integer counts")
    return matrix.astype(np.int64)


def open_data(dataset):
    """Keep original row positions while applying explicit assignment/exclusion inputs."""
    adata, handle = _open_minimal_backed_adata(Path(dataset["h5ad"]))
    try:
        obs = adata.obs.copy()
        obs.index = obs.index.astype(str)
        obs["_row"] = np.arange(len(obs))
        if dataset.get("guide_calls_csv"):
            calls = pd.read_csv(dataset["guide_calls_csv"], dtype="string",
                                usecols=["cell_barcode", "gene_target"])
            require_unique(calls, ["cell_barcode"], "guide assignments")
            require_barcodes(obs.index, calls["cell_barcode"].tolist())
            obs = obs.loc[calls["cell_barcode"]].copy()
            obs[dataset["pert_col"]] = calls["gene_target"].to_numpy()
        settings = dataset.get("statistics", {})
        required = [dataset["pert_col"], settings.get("guide_col", "")]
        required += settings.get("context_cols", [])
        required += [settings[key] for key in ("sample_col", "batch_col") if settings.get(key)]
        if any(not col or col not in obs for col in required):
            raise ValueError(f"Missing declared statistical metadata columns: {required}")
        if obs[required].isna().any().any():
            raise ValueError("Missing labels in statistical metadata; prepare the cohort explicitly")
        if any(obs[col].astype("string").str.strip().eq("").any() for col in required):
            raise ValueError("Empty labels in statistical metadata; prepare the cohort explicitly")
        guide_col = settings["guide_col"]
        obs = obs.loc[~obs[guide_col].astype(str).isin(dataset.get("exclude_guides", []))].copy()
        meta = pd.DataFrame(index=obs.index)
        meta["target"] = obs[dataset["pert_col"]].astype(str)
        meta["guide"] = obs[guide_col].astype(str)
        meta["batch"] = obs[settings["batch_col"]].astype(str) if settings.get("batch_col") else "ALL"
        meta["sample"] = obs[settings["sample_col"]].astype(str) if settings.get("sample_col") else meta["batch"]
        context_cols = settings.get("context_cols", [])
        meta["context"] = (
            obs[context_cols].astype(str).apply(lambda row: json.dumps(row.tolist()), axis=1)
            if context_cols else "ALL"
        )
        meta["_row"] = obs["_row"]
        if meta.empty or dataset["control_label"] not in set(meta["target"]):
            raise ValueError("The declared statistical cohort has no control cells")
        if (meta.groupby("guide", observed=True)["target"].nunique() > 1).any():
            raise ValueError("A guide maps to multiple targets")
        return adata, handle, meta
    except Exception:
        _close_backing_resources(adata, handle)
        raise


def aggregate(adata, handle, membership, read_cells, max_nnz):
    """Stream CSC columns or CSR/dense row blocks; never densify the cell matrix."""
    n_pools, n_cells = membership.shape
    if n_cells != adata.n_obs:
        raise ValueError("Membership does not match the input cell axis")
    encoding = getattr(adata.X, "format", "")
    if encoding == "csc" and handle is not None:
        node = handle["X"]
        data, indices, indptr = [], [], [0]
        for gene in range(adata.n_vars):
            lo, hi = map(int, node["indptr"][gene:gene + 2])
            values = integer_counts(node["data"][lo:hi])
            column = sparse.csc_matrix((values, node["indices"][lo:hi], [0, hi - lo]),
                                       shape=(n_cells, 1))
            counts = (membership @ column).tocsc()
            counts.eliminate_zeros()
            if (counts.data < 0).any():
                raise OverflowError("Pseudobulk integer count overflow")
            data.append(counts.data)
            indices.append(counts.indices)
            indptr.append(indptr[-1] + counts.nnz)
            if indptr[-1] > max_nnz:
                raise ValueError("Pseudobulk nnz exceeds the declared memory budget")
        return sparse.csc_matrix((np.concatenate(data), np.concatenate(indices), indptr),
                                 shape=(n_pools, adata.n_vars))
    counts = sparse.csr_matrix((n_pools, adata.n_vars), dtype=np.int64)
    for start in range(0, n_cells, read_cells):
        stop = min(start + read_cells, n_cells)
        block = integer_counts(adata.X[start:stop, :])
        counts += membership[:, start:stop] @ sparse.csr_matrix(block)
        if counts.nnz > max_nnz:
            raise ValueError("Pseudobulk nnz exceeds the declared memory budget")
        if (counts.data < 0).any():
            raise OverflowError("Pseudobulk integer count overflow")
    return counts.tocsc()


def regular_membership(meta, n_cells):
    keys = ["context", "target", "sample", "batch", "guide"]
    grouped = meta.groupby(keys, sort=True, observed=True)
    codes = grouped.ngroup().to_numpy()
    pools = grouped.size().rename("n_cells").reset_index()
    pools["contrast_id"] = ""
    pools["draw"] = -1
    membership = sparse.csr_matrix(
        (np.ones(len(meta), dtype=np.int64), (codes, meta["_row"])),
        shape=(len(pools), n_cells),
    )
    return membership, pools


def null_guide_sets(frame, control, seed, n_groups):
    guides = sorted(frame.loc[frame["target"].eq(control), "guide"].unique())
    if n_groups == 0:
        return []
    if len(guides) < 3 or math.comb(len(guides), 2) < n_groups:
        raise ValueError("NTC calibration requires enough distinct control-guide pairs; "
                         "provide more controls or explicitly set ntc_null_groups=0")
    rng, selected = np.random.default_rng(seed), set()
    while len(selected) < n_groups:
        selected.add(tuple(sorted(rng.choice(guides, 2, replace=False))))
    return sorted(selected)


def resample_membership(meta, n_cells, control, seed, n_resamples, fraction, null_groups=0):
    """KOLF construction: independent, overlapping draws; controls belong to a contrast."""
    pools, row_indices, cell_indices = [], [], []
    for stratum, frame in meta.groupby(["context", "sample", "batch"], sort=True):
        controls = frame.loc[frame["target"].eq(control), "_row"].to_numpy()
        if len(controls) == 0:
            raise ValueError(f"No same-stratum controls for {stratum}")
        tasks = [(target, cells["_row"].to_numpy(), controls, "target")
                 for target, cells in frame.loc[frame["target"].ne(control)].groupby("target", sort=True)]
        for index, guides in enumerate(null_guide_sets(
                frame, control, stable_seed(seed, *stratum, "null"), null_groups)):
            held_out = frame["target"].eq(control) & frame["guide"].isin(guides)
            remaining = frame["target"].eq(control) & ~frame["guide"].isin(guides)
            tasks.append((f"NTC-null-{index}", frame.loc[held_out, "_row"].to_numpy(),
                          frame.loc[remaining, "_row"].to_numpy(), "ntc_null"))
        for target, positions, comparison_controls, level in tasks:
            size = int(len(positions) * fraction)
            if size > len(comparison_controls):
                raise ValueError(f"Insufficient matched controls for resampling {target}, {stratum}")
            target_seed = stable_seed(seed, *stratum, target)
            rng = np.random.default_rng(target_seed)
            contrast_id = json.dumps([*stratum, target])
            for draw in range(n_resamples):
                for label, candidates in ((target, positions), (control, comparison_controls)):
                    chosen = rng.choice(candidates, size=size, replace=False)
                    pool_index = len(pools)
                    pools.append(dict(context=stratum[0], sample=stratum[1], batch=stratum[2],
                                      target=label, guide="", n_cells=size, contrast_id=contrast_id,
                                      draw=draw, seed=str(target_seed), analysis_level=level))
                    row_indices.extend([pool_index] * size)
                    cell_indices.extend(chosen.tolist())
    if not pools:
        raise ValueError("No perturbations available for KOLF resampling")
    membership = sparse.csr_matrix(
        (np.ones(len(cell_indices), dtype=np.int64), (row_indices, cell_indices)),
        shape=(len(pools), n_cells),
    )
    return membership, pd.DataFrame(pools)


def close_data(adata, handle):
    _close_backing_resources(adata, handle)
